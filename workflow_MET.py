import awkward as ak
import os
import cachetools

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.jets import jet_correction, met_correction_after_jec

from pocket_coffea.lib.deltaR_matching import object_matching, deltaR_matching_nonunique
from custom_cut_functions import *
from custom_functions import *

from params.binning import *
from workflow import QCDBaseProcessor


class METProcessor(QCDBaseProcessor):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

    def apply_object_preselection(self, variation):
        # super().apply_object_preselection(variation)

        print(self.events.Jet.pt_raw, self.events.Jet.pt, self.events.Jet.pt*(1-self.events.Jet.rawFactor))
        cache = cachetools.Cache(np.inf)
        jets_not_calibrated = {}
        jets_calibrated = {}


        if int(os.environ.get("CLOSURE", 0)) == 1:
            for jet_type, jet_coll_name in self.params.jets_calibration.collection[
                self._year
            ].items():
                if "chs" in jet_type or "Puppi" in jet_type:
                    continue
                print(jet_type, jet_coll_name)
                
                # define the pnet reg jet colleciton
                jets_not_calibrated[jet_coll_name] = ak.with_field(
                    self.events["Jet"],
                    ak.where(
                        self.events.Jet.PNetRegPtRawCorr > 0,
                        self.events.Jet.pt_raw
                        * self.events.Jet.PNetRegPtRawCorr
                        * (
                            self.events.Jet.PNetRegPtRawCorrNeutrino
                            if "Neutrino" in jet_coll_name
                            else 1
                        ),
                        self.events.Jet.pt,
                    ),
                    "pt",
                )
                # jets_not_calibrated[jet_coll_name] = ak.with_field(
                #     jets_not_calibrated[jet_coll_name],
                #     jets_not_calibrated[jet_coll_name].pt,
                #     "pt_raw",
                # )
                jets_not_calibrated[jet_coll_name] = ak.with_field(
                    jets_not_calibrated[jet_coll_name],
                    ak.zeros_like(jets_not_calibrated[jet_coll_name].pt),
                    "rawFactor",
                )
                
                jets_not_calibrated[jet_coll_name] = ak.with_field(
                    self.events["Jet"],
                    ak.where(
                        self.events.Jet.PNetRegPtRawCorr > 0,
                        self.events.Jet.mass_raw
                        * self.events.Jet.PNetRegPtRawCorr
                        * (
                            self.events.Jet.PNetRegPtRawCorrNeutrino
                            if "Neutrino" in jet_coll_name
                            else 1
                        ),
                        self.events.Jet.mass,
                    ),
                    "mass",
                )
                # jets_not_calibrated[jet_coll_name] = ak.with_field(
                #     jets_not_calibrated[jet_coll_name],
                #     jets_not_calibrated[jet_coll_name].mass,
                #     "mass_raw",
                # )
                
                jets_calibrated[jet_coll_name] = jet_correction(
                    params=self.params,
                    events=self.events,
                    jets=jets_not_calibrated[jet_coll_name],
                    factory=self.jmefactory,
                    jet_type=jet_type,
                    chunk_metadata={
                        "year": self._year,
                        "isMC": self._isMC,
                        "era": self._era,
                    },
                    cache=cache,
                )
                print(self.events.Jet.PNetRegPtRawCorr, self.events.Jet.PNetRegPtRawCorrNeutrino)
                print(jets_not_calibrated[jet_coll_name].pt, jets_calibrated[jet_coll_name].pt, jets_not_calibrated[jet_coll_name].eta)

                self.events["Jet"] = ak.with_field(
                    self.events["Jet"],
                    jets_calibrated[jet_coll_name].pt,
                    f"pt{jet_coll_name.split('Jet')[1]}",
                )
                self.events["Jet"] = ak.with_field(
                    self.events["Jet"],
                    jets_calibrated[jet_coll_name].mass,
                    f"mass{jet_coll_name.split('Jet')[1]}",
                )
                
                print(getattr(self.events.Jet, f"pt{jet_coll_name.split('Jet')[1]}"))


            # TODO: correct the MET with the new

    def process_extra_after_presel(self, variation) -> ak.Array:
        # super().process_extra_after_presel(variation)
        pass
