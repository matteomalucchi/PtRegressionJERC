import awkward as ak
import os
import cachetools
import numpy as np
import vector

vector.register_awkward()

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.jets import jet_correction, met_correction_after_jec
from pocket_coffea.lib.leptons import lepton_selection, get_dilepton
from pocket_coffea.lib.deltaR_matching import object_matching, deltaR_matching_nonunique

from configs.jme.workflow import QCDBaseProcessor
from configs.jme.custom_cut_functions import jet_selection_nopu
from utils.basic_functions import add_fields


class METProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.only_physical_jet = self.workflow_options["only_physical_jet"]
        self.rescale_MET_with_regressed_pT = self.workflow_options[
            "rescale_MET_with_regressed_pT"
        ]
        self.jec_pt_threshold = self.workflow_options["jec_pt_threshold"]
        self.consider_all_jets = self.workflow_options["consider_all_jets"]
        self.add_corr_t1_met_jets = self.workflow_options["add_corr_t1_met_jets"]
        self.jet_type1_selections = self.workflow_options["jet_type1_selections"]

    def add_GenMET_plus_neutrino(self):
        # Add the neutrinos to the GetJets to compute the MET with neutrinos
        neutrinos = self.events["GenPart"][
            (abs(self.events.GenPart.pdgId) == 12)
            | (abs(self.events.GenPart.pdgId) == 14)
            | (abs(self.events.GenPart.pdgId) == 16)
        ]
        self.events["GenJetPlusNeutrino"] = QCDBaseProcessor.add_neutrinos_to_genjets(
            self, self.events["GenJetGood"], neutrinos
        )

        GenMETPlusNeutrino = met_correction_after_jec(
            self.events,
            "GenMET",
            self.events.GenJetGood,
            self.events.GenJetPlusNeutrino,
        )
        self.events["GenMETPlusNeutrino"] = ak.with_field(
            self.events["GenMET"], GenMETPlusNeutrino["pt"], "pt"
        )
        self.events["GenMETPlusNeutrino"] = ak.with_field(
            self.events["GenMETPlusNeutrino"], GenMETPlusNeutrino["phi"], "phi"
        )

    def process_extra_after_skim(self):
        if "pt_raw" not in self.events["Jet"].fields:
            self.events["Jet"] = ak.with_field(
                self.events["Jet"],
                self.events["Jet"].pt * (1 - self.events["Jet"].rawFactor),
                "pt_raw",
            )
            self.events["Jet"] = ak.with_field(
                self.events["Jet"],
                self.events["Jet"].mass * (1 - self.events["Jet"].rawFactor),
                "mass_raw",
            )
            
        if self.add_corr_t1_met_jets:
            # consider the lowpt_jets collection
            lowpt_jets = self.events["CorrT1METJet"]
            lowpt_jets = ak.with_field(lowpt_jets, lowpt_jets.rawPt, "pt")
            lowpt_jets = ak.with_field(lowpt_jets, lowpt_jets.rawPt, "pt_raw")
            lowpt_jets = ak.with_field(
                lowpt_jets, ak.zeros_like(lowpt_jets.pt, dtype=np.float32), "mass"
            )
            lowpt_jets = ak.with_field(
                lowpt_jets, ak.zeros_like(lowpt_jets.pt, dtype=np.float32), "mass_raw"
            )
            lowpt_jets = ak.with_field(
                lowpt_jets, ak.zeros_like(lowpt_jets.pt, dtype=np.float32), "rawFactor"
            )
            lowpt_jets = ak.with_field(lowpt_jets, ak.ones_like(lowpt_jets.pt, dtype=np.float32), "PNetRegPtRawCorr")
            lowpt_jets = ak.with_field(lowpt_jets, ak.ones_like(lowpt_jets.pt, dtype=np.float32), "PNetRegPtRawCorrNeutrino")
            lowpt_jets = ak.with_field(lowpt_jets, ak.ones_like(lowpt_jets.pt, dtype=np.float32), "btagPNetB")
            lowpt_jets = ak.with_field(lowpt_jets, ak.ones_like(lowpt_jets.pt, dtype=np.float32), "btagPNetCvL")
            
            lowpt_jets=add_fields(lowpt_jets, "all", four_vec="Momentum4D")

            # concatenate the two collections
            self.events["JetGood"] = ak.concatenate(
                [self.events["Jet"], lowpt_jets], axis=1
            )
        elif self.consider_all_jets:
            self.events["JetGood"] = ak.copy(self.events["Jet"])
        else:
            # keep only jets with pt_raw > 15 GeV and |eta| < 4.7
            # TODO  add here the cut on pt only for the regression
            self.events["JetGood"] = jet_selection_nopu(
                self.events, "Jet", self.params, "pt_raw"
            )

            if self.only_physical_jet:
                physisical_jet_mask = (
                    self.events["JetGood"].pt_raw * np.cosh(self.events["JetGood"].eta)
                    < (13.6 * 1000) / 2
                )
                self.events["JetGood"] = self.events["JetGood"][physisical_jet_mask]

            # consider only jets with regressed pt > 0
            reg_mask= self.events["JetGood"].PNetRegPtRawCorr>0
            self.events["JetGood"] = self.events["JetGood"][reg_mask]
            
        # Create extra Jet collections for calibration
        self.events["JetGoodJEC"] = ak.copy(self.events["JetGood"])
        self.events["JetGoodPNet"] = ak.copy(self.events["JetGood"])
        self.events["JetGoodPNetPlusNeutrino"] = ak.copy(self.events["JetGood"])

    def apply_object_preselection(self, variation):
        self.events["GenJetGood"] = self.events.GenJet[
            self.events.GenJet.pt > self.params.object_preselection["GenJet"]["pt"]
        ]

        jets_raw = ak.zip(
            {
                "pt": self.events.JetGood.pt_raw,
                "eta": self.events.JetGood.eta,
                "phi": self.events.JetGood.phi,
                "mass": self.events.JetGood.mass_raw,
            },
            with_name="Momentum4D",
        )

        self.met_branches = ["RawPuppiMET", "PuppiMET"]

        for jet_coll_name in [
            "JetGood",
            "JetGoodJEC",
            "JetGoodPNet",
            "JetGoodPNetPlusNeutrino",
        ]:

            self.events[jet_coll_name]=add_fields(self.events[jet_coll_name], four_vec="Momentum4D")
            if "PNet" in jet_coll_name:
                # Correct MET with MC Truth only jets with pt reg > 15
                corr_reg_pt_mask = (
                    self.events[jet_coll_name].pt_raw > self.jec_pt_threshold
                )
                # compute pt
                self.events[jet_coll_name] = ak.with_field(
                    self.events[jet_coll_name],
                    ak.where(
                        corr_reg_pt_mask,
                        self.events[jet_coll_name].pt,
                        self.events[jet_coll_name].pt_raw,
                    ),
                    "pt",
                )
                # compute mass
                self.events[jet_coll_name] = ak.with_field(
                    self.events[jet_coll_name],
                    ak.where(
                        corr_reg_pt_mask,
                        self.events[jet_coll_name].mass,
                        self.events[jet_coll_name].mass_raw,
                    ),
                    "mass",
                )
                # compute px and py
                self.events[jet_coll_name] = ak.with_field(
                    self.events[jet_coll_name],
                    self.events[jet_coll_name].pt
                    * np.cos(self.events[jet_coll_name].phi),
                    "px",
                )
                self.events[jet_coll_name] = ak.with_field(
                    self.events[jet_coll_name],
                    self.events[jet_coll_name].pt
                    * np.sin(self.events[jet_coll_name].phi),
                    "py",
                )

            if self.rescale_MET_with_regressed_pT:
                for met_branch, jet_coll in zip(
                    ["RawPuppiMET", "PuppiMET"],
                    [jets_raw, self.events.JetGood],
                ):

                    new_MET = met_correction_after_jec(
                        self.events,
                        met_branch,
                        jet_coll,
                        self.events[jet_coll_name],
                    )
                    jet_coll_suffix = jet_coll_name.split("JetGood")[-1]
                    new_met_branch = f"{met_branch}-Type1{jet_coll_suffix}"

                    self.events[new_met_branch] = ak.zip(
                        {
                            "pt": new_MET["pt"],
                            "phi": new_MET["phi"],
                        },
                    )
                    self.met_branches.append(new_met_branch)

        self.events["MuonGood"] = lepton_selection(self.events, "Muon", self.params)
        self.events["ElectronGood"] = lepton_selection(
            self.events, "Electron", self.params
        )
        self.events["ll"] = get_dilepton(None, self.events.MuonGood)

    def get_hadronic_recoil(self, lepton_coll, MET_coll):
        met = self.events[MET_coll]
        leptons_px = ak.sum(
            self.events[lepton_coll].pt * np.cos(self.events[lepton_coll].phi), axis=1
        )
        leptons_py = ak.sum(
            self.events[lepton_coll].pt * np.sin(self.events[lepton_coll].phi), axis=1
        )
        hadronic_recoil_px = -(met.pt * np.cos(met.phi) + leptons_px)
        hadronic_recoil_py = -(met.pt * np.sin(met.phi) + leptons_py)

        hadronic_recoil_pt = np.sqrt(hadronic_recoil_px**2 + hadronic_recoil_py**2)
        hadronic_recoil_phi = np.arctan2(hadronic_recoil_py, hadronic_recoil_px)

        hadronic_recoil = ak.zip(
            {
                "px": hadronic_recoil_px,
                "py": hadronic_recoil_py,
                "pt": hadronic_recoil_pt,
                "phi": hadronic_recoil_phi,
            },
            with_name="Momentum2D",
        )

        return hadronic_recoil

    def compute_projections_qT(self, hadronic_recoil_coll):
        v_qT = self.events["qT"]

        # hadronic recoil
        u_vec = self.events[hadronic_recoil_coll]
        # compute dot product of v_qT and MET
        response = v_qT.dot(-u_vec) / v_qT.dot(v_qT)

        v_paral_predict = v_qT.dot(u_vec) / v_qT.rho

        # subtract the module of qT because the projection
        # of u is in the direction opposite of qT
        u_paral_predict = v_paral_predict + v_qT.rho

        # vector 90 degrees to the left of v_qT
        v_qT_perp = ak.zip(
            {
                "x": -v_qT.y,
                "y": v_qT.x,
            },
            with_name="Vector2D",
        )

        u_perp_predict = v_qT_perp.dot(u_vec) / v_qT_perp.rho

        return u_perp_predict, u_paral_predict, response

    def process_extra_after_presel(self, variation) -> ak.Array:
        self.events["qT"] = ak.zip(
            {"x": self.events["ll"].px, "y": self.events["ll"].py},
            with_name="Vector2D",
        )
        # substract leptons from MET
        for MET_coll in self.met_branches:
            lepton_coll = "MuonGood"
            self.events[f"u{MET_coll}"] = self.get_hadronic_recoil(
                lepton_coll, MET_coll
            )

            u_perp_predict, u_paral_predict, response = self.compute_projections_qT(
                f"u{MET_coll}"
            )
            self.events[f"u{MET_coll}"] = ak.with_field(
                self.events[f"u{MET_coll}"],
                u_perp_predict,
                "u_perp_predict",
            )
            self.events[f"u{MET_coll}"] = ak.with_field(
                self.events[f"u{MET_coll}"],
                u_paral_predict,
                "u_paral_predict",
            )
            self.events[f"u{MET_coll}"] = ak.with_field(
                self.events[f"u{MET_coll}"], response, "response"
            )

    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nJetGood"] = ak.num(self.events.JetGood)
