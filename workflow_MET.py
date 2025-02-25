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
from workflow  import QCDBaseProcessor


class METProcessor(QCDBaseProcessor):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)



    def apply_object_preselection(self, variation):
        super().apply_object_preselection(variation)
        cache = cachetools.Cache(np.inf)
        for jet_type, jet_coll_name in self.params.jets_calibration.collection[self._year].items():
            
            # NOTE: consider the jetgood when and not all jets
            self.events["JetPNet"] = jet_correction(
                params=self.params,
                events=self.events,
                jets=self.events[jet_coll_name],
                factory=self.jmefactory,
                jet_type = jet_type,
                chunk_metadata={
                    "year": self._year,
                    "isMC": self._isMC,
                    "era": self._era,
                },
                cache=cache
            )
            
            # TODO: correct the MET with the new

    def process_extra_after_presel(self, variation) -> ak.Array:
        super().process_extra_after_presel(variation)
