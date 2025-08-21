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


class METProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.only_physical_jet = self.workflow_options["only_physical_jet"]
        self.rescale_MET_with_regressed_pT = self.workflow_options[
            "rescale_MET_with_regressed_pT"
        ]
        self.jec_pt_threshold = self.workflow_options["jec_pt_threshold"]

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
        
        # Create extra Jet collections for calibration
        self.events["JetGoodJEC"] = ak.copy(self.events["JetGood"])
        self.events["JetGoodPtReg"] = ak.copy(self.events["JetGood"])
        self.events["JetGoodPtRegPlusNeutrino"] = ak.copy(self.events["JetGood"])

    def apply_object_preselection(self, variation):
        # if "pt_raw" not in self.events["Jet"].fields:
        #     self.events["Jet"] = ak.with_field(
        #         self.events["Jet"],
        #         self.events["Jet"].pt * (1 - self.events["Jet"].rawFactor),
        #         "pt_raw",
        #     )
        #     self.events["Jet"] = ak.with_field(
        #         self.events["Jet"],
        #         self.events["Jet"].mass * (1 - self.events["Jet"].rawFactor),
        #         "mass_raw",
        #     )

        # # keep only jets with pt_raw > 15 GeV and |eta| < 4.7
        # for suffix in ["", "JEC", "PtReg", "PtRegPlusNeutrino"]:
        #     self.events[f"JetGood{suffix}"] = jet_selection_nopu(
        #         self.events, f"Jet{suffix}", self.params, "pt_raw"
        #     )
        #     if self.only_physical_jet:
        #         physisical_jet_mask = (
        #             self.events[f"JetGood{suffix}"].pt_raw * np.cosh(self.events[f"JetGood{suffix}"].eta)
        #             < (13.6 * 1000) / 2
        #         )
        #         self.events[f"JetGood{suffix}"] = self.events[f"JetGood{suffix}"][physisical_jet_mask]
        

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

        # cache = cachetools.Cache(np.inf)
        # jets_dict = {}
        # jet_calib_params = self.params.jets_calibration

        for jet_coll_name in ["JetGood", "JetGoodJEC", "JetGoodPtReg", "JetGoodPtRegPlusNeutrino"]:
            # for _, jet_coll_name in jet_calib_params.collection[self._year].items():
            # if "chs" in jet_type or "Puppi" in jet_type:
            #     continue

            # if "PNet" in jet_coll_name:
            #     # For regression
            #     # define the pnet reg jet colleciton
            #     jets_dict[jet_coll_name] = ak.with_field(
            #         self.events["JetGood"],
            #         self.events.JetGood.pt_raw
            #         * self.events.JetGood.PNetRegPtRawCorr
            #         * (
            #             self.events.JetGood.PNetRegPtRawCorrNeutrino
            #             if "Neutrino" in jet_coll_name
            #             else 1
            #         ),
            #         "pt",
            #     )
            #     jets_dict[jet_coll_name] = ak.with_field(
            #         jets_dict[jet_coll_name],
            #         jets_dict[jet_coll_name].pt,
            #         "pt_raw",
            #     )
            #     jets_dict[jet_coll_name] = ak.with_field(
            #         jets_dict[jet_coll_name],
            #         ak.zeros_like(jets_dict[jet_coll_name].pt),
            #         "rawFactor",
            #     )

            #     jets_dict[jet_coll_name] = ak.with_field(
            #         jets_dict[jet_coll_name],
            #         jets_dict[jet_coll_name].mass_raw
            #         * jets_dict[jet_coll_name].PNetRegPtRawCorr
            #         * (
            #             jets_dict[jet_coll_name].PNetRegPtRawCorrNeutrino
            #             if "Neutrino" in jet_coll_name
            #             else 1
            #         ),
            #         "mass",
            #     )
            #     jets_dict[jet_coll_name] = ak.with_field(
            #         jets_dict[jet_coll_name],
            #         jets_dict[jet_coll_name].mass,
            #         "mass_raw",
            #     )
            # else:
            #     # For nominal jets
            #     jets_dict[jet_coll_name] = ak.copy(self.events["JetGood"])

            # # Calibrate the jets
            # jets_calib_dict = {}
            # # always apply the JEC to the regressed pt
            # if True or jet_calib_params.apply_jec_nominal[self._year]:

            #     jets_calib_dict[jet_coll_name] = jet_correction(
            #         params=self.params,
            #         events=self.events,
            #         jets=jets_dict[jet_coll_name],
            #         factory=self.jmefactory,
            #         jet_type=jet_type,
            #         chunk_metadata={
            #             "year": self._year,
            #             "isMC": self._isMC,
            #             "era": self._era,
            #         },
            #         cache=cache,
            #     )
            # else:
            #     jets_calib_dict[jet_coll_name] = jets_dict[jet_coll_name]

            # jet_coll_suffix = jet_coll_name.split("Jet")[-1]

            # Correct MET with MC Truth only jets with pt reg > 15
            # corr_reg_pt_mask = (
            #     jets_calib_dict[jet_coll_name].pt_raw > self.jec_pt_threshold
            # )
            # # compute pt
            # self.events[f"JetGood{jet_coll_suffix}"] = ak.with_field(
            #     jets_calib_dict[jet_coll_name],
            #     ak.where(
            #         corr_reg_pt_mask,
            #         jets_calib_dict[jet_coll_name].pt,
            #         jets_calib_dict[jet_coll_name].pt_raw,
            #     ),
            #     "pt",
            # )
            
            
            
            # Correct MET with MC Truth only jets with pt reg > 15
            corr_reg_pt_mask = (
                self.events[jet_coll_name].pt_raw > self.jec_pt_threshold
            )
            breakpoint()
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
            breakpoint()

            # # compute px and py
            # self.events[f"JetGood{jet_coll_suffix}"] = ak.with_field(
            #     self.events[f"JetGood{jet_coll_suffix}"],
            #     self.events[f"JetGood{jet_coll_suffix}"].pt
            #     * np.cos(self.events[f"JetGood{jet_coll_suffix}"].phi),
            #     "px",
            # )
            # self.events[f"JetGood{jet_coll_suffix}"] = ak.with_field(
            #     self.events[f"JetGood{jet_coll_suffix}"],
            #     self.events[f"JetGood{jet_coll_suffix}"].pt
            #     * np.sin(self.events[f"JetGood{jet_coll_suffix}"].phi),
            #     "py",
            # )


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
                    jet_coll_suffix = jet_coll_name.split("Jet")[-1]

                    self.events[f"{met_branch}-{jet_coll_suffix}"] = ak.zip(
                        {
                            "pt": new_MET["pt"],
                            "phi": new_MET["phi"],
                        },
                    )

        self.events["MuonGood"] = lepton_selection(self.events, "Muon", self.params)
        self.events["ElectronGood"] = lepton_selection(
            self.events, "Electron", self.params
        )
        self.events["ll"] = get_dilepton(None, self.events.MuonGood)

    def subtract_leptons_from_MET(self, lepton_coll, MET_coll):
        met = self.events[MET_coll]
        leptons_px = ak.sum(
            self.events[lepton_coll].pt * np.cos(self.events[lepton_coll].phi), axis=1
        )
        leptons_py = ak.sum(
            self.events[lepton_coll].pt * np.sin(self.events[lepton_coll].phi), axis=1
        )
        new_met_px = met.pt * np.cos(met.phi) + leptons_px
        new_met_py = met.pt * np.sin(met.phi) + leptons_py

        new_met_pt = np.sqrt(new_met_px**2 + new_met_py**2)
        new_met_phi = np.arctan2(new_met_py, new_met_px)

        return ak.zip(
            {"px": new_met_px, "py": new_met_py, "pt": new_met_pt, "phi": new_met_phi},
            with_name="Momentum2D",
        )

    def compute_projections_qT(self, MET_coll):
        v_qT = self.events["qT"]

        # hadronic recoil
        u_vec = -self.events[MET_coll]
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
        for MET_coll in [
            "RawPuppiMET",
            "RawPuppiMET-JEC",
            "RawPuppiMET-PNet",
            "RawPuppiMET-PNetPlusNeutrino",
            "PuppiMET",
            "PuppiMET-JEC",
            "PuppiMET-PNet",
            "PuppiMET-PNetPlusNeutrino",
            # "GenMET",
            # "GenMETPlusNeutrino",
        ]:
            lepton_coll = "MuonGood"
            self.events[f"{MET_coll}_{lepton_coll}"] = self.subtract_leptons_from_MET(
                lepton_coll, MET_coll
            )

            u_perp_predict, u_paral_predict, response = self.compute_projections_qT(
                f"{MET_coll}_{lepton_coll}"
            )
            self.events[f"{MET_coll}_{lepton_coll}"] = ak.with_field(
                self.events[f"{MET_coll}_{lepton_coll}"],
                u_perp_predict,
                "u_perp_predict",
            )
            self.events[f"{MET_coll}_{lepton_coll}"] = ak.with_field(
                self.events[f"{MET_coll}_{lepton_coll}"],
                u_paral_predict,
                "u_paral_predict",
            )
            self.events[f"{MET_coll}_{lepton_coll}"] = ak.with_field(
                self.events[f"{MET_coll}_{lepton_coll}"], response, "response"
            )

    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nJetGood"] = ak.num(self.events.JetGood)
