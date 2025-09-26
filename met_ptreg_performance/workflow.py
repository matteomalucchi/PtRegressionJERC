import awkward as ak
import os
import cachetools
import numpy as np
import vector

vector.register_awkward()

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.jets import (
    jet_correction,
    met_correction_after_jec,
    jet_type1_selection,
)
from pocket_coffea.lib.leptons import lepton_selection, get_dilepton
from pocket_coffea.lib.deltaR_matching import object_matching, deltaR_matching_nonunique

from configs.jme.workflow import QCDBaseProcessor
from configs.jme.custom_cut_functions import jet_selection_nopu
from utils.basic_functions import add_fields
from configs.MET_studies.custom_cuts_functions import muon_selection_custom


class METProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.only_physical_jet = self.workflow_options["only_physical_jet"]
        self.rescale_MET_with_regressed_pT = self.workflow_options[
            "rescale_MET_with_regressed_pT"
        ]
        self.jec_pt_threshold = self.workflow_options["jec_pt_threshold"]
        self.consider_all_jets = self.workflow_options["consider_all_jets"]

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

        # compute EmEF for the jets
        self.events["Jet"] = ak.with_field(
            self.events["Jet"],
            self.events["Jet"].chEmEF + self.events["Jet"].neEmEF,
            "EmEF",
        )
        if self.consider_all_jets:
            self.events["JetGood"] = ak.copy(self.events["Jet"])
        else:
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
            # keep only jets with pt_raw > 15 GeV and |eta| < 4.7
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
            reg_mask = self.events["JetGood"].PNetRegPtRawCorr > 0
            self.events["JetGood"] = self.events["JetGood"][reg_mask]

        # Create extra Jet collections for calibration
        self.events["JetGoodJEC"] = ak.copy(self.events["JetGood"])
        self.events["JetGoodPNet"] = ak.copy(self.events["JetGood"])
        self.events["JetGoodPNetPlusNeutrino"] = ak.copy(self.events["JetGood"])

        # Jets for type 1 met correction
    
        # consider the lowpt jets collection for type 1 met correction
        self.events["JetLowPtMuonSubtr"] = ak.copy(self.events["CorrT1METJet"])
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"], self.events["JetLowPtMuonSubtr"].rawPt, "pt"
        )
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"],
            ak.zeros_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32),
            "mass",
        )
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"],
            ak.zeros_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32),
            "rawFactor",
        )
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"],
            ak.ones_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32),
            "PNetRegPtRawCorr",
        )
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"],
            ak.ones_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32),
            "PNetRegPtRawCorrNeutrino",
        )
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"],
            ak.zeros_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32) - 1,
            "btagPNetB",
        )
        self.events["JetLowPtMuonSubtr"] = ak.with_field(
            self.events["JetLowPtMuonSubtr"],
            ak.zeros_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32) - 1,
            "btagPNetCvL",
        )
        if "EmEF" not in self.events["JetLowPtMuonSubtr"].fields:
            # EmEF is needed for the jet cleaning in the type1 met correction
            self.events["JetLowPtMuonSubtr"] = ak.with_field(
                self.events["JetLowPtMuonSubtr"],
                ak.zeros_like(self.events["JetLowPtMuonSubtr"].pt, dtype=np.float32) - 1,
                "EmEF",
            )
        self.events["JetLowPtMuonSubtr"] = add_fields(
            self.events["JetLowPtMuonSubtr"], "all", four_vec="Momentum4D"
        )

        self.events["JetMuonSubtr"] = ak.copy(self.events["Jet"])
        
        for jet_coll in ["JetMuonSubtr", "JetLowPtMuonSubtr"]:
            # Change the rawFactor of the JetMuonSubtr and JetLowPtMuonSubtr collection to account for muon subtraction
            # This is not exactly the correct thing to do. We shouldn't be using the full
            # jet eta and phi but the eta and phi of muon-less jet p4. Not possible with
            # current NanoAODs.
            # This is the factor that is used to get the raw muon-less jet pT and mass
            self.events[jet_coll] = ak.with_field(
                self.events[jet_coll],
                1
                - (
                    (1 - self.events[jet_coll].rawFactor)
                    * (1 - self.events[jet_coll].muonSubtrFactor)
                ),
                "rawFactor",
            )

            if "muonSubtrDeltaEta" in self.events[jet_coll].fields:
                self.events[jet_coll] = ak.with_field(
                    self.events[jet_coll],
                    self.events[jet_coll].muonSubtrDeltaEta
                    + self.events[jet_coll].eta,
                    "eta",
                )
            if "muonSubtrDeltaPhi" in self.events[jet_coll].fields:
                self.events[jet_coll] = ak.with_field(
                    self.events[jet_coll],
                    self.events[jet_coll].muonSubtrDeltaPhi
                    + self.events[jet_coll].phi,
                    "phi",
                )
                
    def apply_object_preselection(self, variation):
        self.events["GenJetGood"] = self.events.GenJet[
            self.events.GenJet.pt > self.params.object_preselection["GenJet"]["pt"]
        ]

        jets_raw_coll = ak.zip(
            {
                "pt": self.events.JetGood.pt_raw,
                "eta": self.events.JetGood.eta,
                "phi": self.events.JetGood.phi,
                "mass": self.events.JetGood.mass_raw,
            },
            with_name="Momentum4D",
        )
        
        self.events["JetCorrMET"] = ak.concatenate(
            [ak.copy(self.events["JetMuonSubtr"]), self.events["JetLowPtMuonSubtr"]],
            axis=1,
        )

        self.events["JetGoodCorrMET"] = jet_type1_selection(
            self.events, "JetCorrMET", self.params
        )

        # for the JetGoodCorrMET we use the raw mass
        jets_corrmet_raw_coll = ak.zip(
            {
                "pt": self.events.JetGoodCorrMET.pt_raw,
                "eta": self.events.JetGoodCorrMET.eta,
                "phi": self.events.JetGoodCorrMET.phi,
                "mass": self.events.JetGoodCorrMET.mass_raw,
            },
            with_name="Momentum4D",
        )

        self.met_branches = ["RawPuppiMET", "PuppiMET"]

        for jet_coll_name, jets_raw in zip(
            [
                "JetGood",
                "JetGoodJEC",
                "JetGoodCorrMET",
                "JetGoodPNet",
                "JetGoodPNetPlusNeutrino",
            ],
            [
                jets_raw_coll,
                jets_raw_coll,
                jets_corrmet_raw_coll,
                jets_raw_coll,
                jets_raw_coll,
            ],
        ):

            self.events[jet_coll_name] = add_fields(
                self.events[jet_coll_name], four_vec="Momentum4D"
            )
            if "PNet" in jet_coll_name:
                continue

                # TODO: this is probably wrong for the regression
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

        self.events["MuonGood"] = muon_selection_custom(self.events, self.params)
        self.events["ElectronGood"] = lepton_selection(
            self.events, "Electron", self.params
        )

        # di lepton is composed by the 2 leading muons
        self.events["ll"] = get_dilepton(None, self.events.MuonGood)

    def get_hadronic_recoil(self, qT, MET_coll):
        met = self.events[MET_coll]

        hadronic_recoil_px = -(met.pt * np.cos(met.phi) + qT.x)
        hadronic_recoil_py = -(met.pt * np.sin(met.phi) + qT.y)

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
        for MET_coll in self.met_branches:
            # substract leptons from MET
            self.events[f"u{MET_coll}"] = self.get_hadronic_recoil(
                self.events["qT"], MET_coll
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
