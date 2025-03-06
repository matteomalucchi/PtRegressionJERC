import awkward as ak
import os
import cachetools
import numpy as np

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.jets import jet_correction, met_correction_after_jec
from pocket_coffea.lib.deltaR_matching import object_matching, deltaR_matching_nonunique

from configs.jme.workflow import QCDBaseProcessor
from configs.jme.custom_cut_functions import jet_selection_nopu


class METProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.only_physisical_jet = self.workflow_options["only_physisical_jet"]


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
        print("GenMet", self.events["GenMET"].pt)
        print("GenMetPlusNeutrino", self.events["GenMETPlusNeutrino"].pt)

    def apply_object_preselection(self, variation):
        # super().apply_object_preselection(variation)
        self.events["JetPuppiMET"] = self.events["Jet"]
        self.events["JetPuppiMET"] = jet_selection_nopu(
            self.events, "JetPuppiMET", self.params, "pt_raw"
        )
        self.events["GenJetGood"] =  self.events.GenJet[self.events.GenJet.pt>self.params.object_preselection["GenJet"]["pt"]]
        

        if self.only_physisical_jet:
            physisical_jet_mask = (
                self.events.JetPuppiMET.pt_raw * np.cosh(self.events.JetPuppiMET.eta)
                < (13.6 * 1000) / 2
            )
            self.events["JetPuppiMET"] = self.events.JetPuppiMET[physisical_jet_mask]

        cache = cachetools.Cache(np.inf)
        jets_dict = {}

        jet_calib_params = self.params.jets_calibration
        for jet_type, jet_coll_name in jet_calib_params.collection[self._year].items():
            if "chs" in jet_type or "Puppi" in jet_type:
                continue
            
            
            # define the pnet reg jet colleciton
            jets_dict[jet_coll_name] = ak.with_field(
                self.events["JetPuppiMET"],
                ak.where(
                    self.events.JetPuppiMET.PNetRegPtRawCorr > 0,
                    self.events.JetPuppiMET.pt_raw
                    * self.events.JetPuppiMET.PNetRegPtRawCorr
                    * (
                        self.events.JetPuppiMET.PNetRegPtRawCorrNeutrino
                        if "Neutrino" in jet_coll_name
                        else 1
                    ),
                    self.events.JetPuppiMET.pt,
                ),
                "pt",
            )
            jets_dict[jet_coll_name] = ak.with_field(
                jets_dict[jet_coll_name],
                jets_dict[jet_coll_name].pt,
                "pt_raw",
            )
            jets_dict[jet_coll_name] = ak.with_field(
                jets_dict[jet_coll_name],
                ak.zeros_like(jets_dict[jet_coll_name].pt),
                "rawFactor",
            )
            
            jets_dict[jet_coll_name] = ak.with_field(
                jets_dict[jet_coll_name],
                ak.where(
                    jets_dict[jet_coll_name].PNetRegPtRawCorr > 0,
                    jets_dict[jet_coll_name].mass_raw
                    * jets_dict[jet_coll_name].PNetRegPtRawCorr
                    * (
                        jets_dict[jet_coll_name].PNetRegPtRawCorrNeutrino
                        if "Neutrino" in jet_coll_name
                        else 1
                    ),
                    jets_dict[jet_coll_name].mass,
                ),
                "mass",
            )
            
            #Calibrate the jets
            jets_calib_dict={}
            if (
                jet_calib_params.apply_jec_nominal[self._year]
            ):
                print(jet_type, jet_coll_name)
                print(
                    "PNETReg before correction",
                    jets_dict[jet_coll_name].pt,
                    jets_dict[jet_coll_name].mass,
                    jets_dict[jet_coll_name].eta,
                )

                jets_calib_dict[jet_coll_name] = jet_correction(
                    params=self.params,
                    events=self.events,
                    jets=jets_dict[jet_coll_name],
                    factory=self.jmefactory,
                    jet_type=jet_type,
                    chunk_metadata={
                        "year": self._year,
                        "isMC": self._isMC,
                        "era": self._era,
                    },
                    cache=cache,
                )
                print(
                    "PNETReg after correction",
                    jets_calib_dict[jet_coll_name].pt,
                    jets_calib_dict[jet_coll_name].mass,
                    jets_calib_dict[jet_coll_name].eta,
                )
            else:
                jets_calib_dict[jet_coll_name] = jets_dict[jet_coll_name]

            jet_coll_suffix = jet_coll_name.split("Jet")[-1]
            print("jet_coll_suffix", jet_coll_suffix)
         
            # compute pt
            reg_pt_mask=jets_dict[jet_coll_name].pt > 15
            self.events[f"JetPuppiMET{jet_coll_suffix}"] = ak.with_field(
                jets_calib_dict[jet_coll_name],
                ak.where(
                    reg_pt_mask,
                    jets_calib_dict[jet_coll_name].pt,
                    jets_dict[jet_coll_name].pt,
                ),
                "pt"
            )
            # compute px and py
            self.events[f"JetPuppiMET{jet_coll_suffix}"] = ak.with_field(
                self.events[f"JetPuppiMET{jet_coll_suffix}"],
                self.events[f"JetPuppiMET{jet_coll_suffix}"].pt * np.cos(self.events[f"JetPuppiMET{jet_coll_suffix}"].phi),
                "px",
            )
            self.events[f"JetPuppiMET{jet_coll_suffix}"] = ak.with_field(
                self.events[f"JetPuppiMET{jet_coll_suffix}"],
                self.events[f"JetPuppiMET{jet_coll_suffix}"].pt * np.sin(self.events[f"JetPuppiMET{jet_coll_suffix}"].phi),
                "py",
            )
            # compute mass
            self.events[f"JetPuppiMET{jet_coll_suffix}"] = ak.with_field(
                self.events[f"JetPuppiMET{jet_coll_suffix}"],
                ak.where(
                    reg_pt_mask,
                    jets_calib_dict[jet_coll_name].mass,
                    jets_dict[jet_coll_name].mass,
                ),
                "mass",
            )

            print("jetpuppimet pt", self.events.JetPuppiMET.pt)
            print("jetpuppimetPNet pt", self.events[f"JetPuppiMET{jet_coll_suffix}"].pt)
            
            print("jetpuppimet mass", self.events.JetPuppiMET.mass)
            print("jetpuppimetPNet mass", self.events[f"JetPuppiMET{jet_coll_suffix}"].mass)

            if jet_calib_params.rescale_MET[self._year]:
                met_branch = jet_calib_params.rescale_MET_branch[self._year]
                print(
                    met_branch, self.events[met_branch].pt, self.events[met_branch].phi
                )
                print(
                    "JetPuppiMET px",
                    self.events["JetPuppiMET"].px,
                    "py",
                    self.events["JetPuppiMET"].py,
                )
                print(
                    "jet_dict px",
                    self.events[f"JetPuppiMET{jet_coll_suffix}"].px,
                    "py",
                    self.events[f"JetPuppiMET{jet_coll_suffix}"].py,
                )

                new_MET = met_correction_after_jec(
                    self.events,
                    met_branch,
		            self.events.JetPuppiMET,
                    self.events[f"JetPuppiMET{jet_coll_suffix}"],
                )
                self.events[f"{met_branch}{jet_coll_suffix}"] = ak.with_field(
                    self.events[met_branch], new_MET["pt"], "pt"
                )
                self.events[f"{met_branch}{jet_coll_suffix}"] = ak.with_field(
                    self.events[f"{met_branch}{jet_coll_suffix}"], new_MET["phi"], "phi"
                )
                print(
                    f"{met_branch}{jet_coll_suffix}",
                    self.events[f"{met_branch}{jet_coll_suffix}"].pt,
                    self.events[f"{met_branch}{jet_coll_suffix}"].phi,
                )


        self.add_GenMET_plus_neutrino()
        
    def process_extra_after_presel(self, variation) -> ak.Array:
        # super().process_extra_after_presel(variation)
        pass
    
    def count_objects(self, variation):
        pass