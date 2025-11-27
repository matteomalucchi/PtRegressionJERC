import awkward as ak
import copy
from pocket_coffea.lib.cut_functions import get_JetVetoMap_Mask


def dimuon(events, params, year, sample, **kwargs):

    # Masks for same-flavor (SF) and opposite-sign (OS)
    SF = (events.nMuonGood >= 2) & (events.nElectronGood == 0)
    OS = events.ll.charge == 0

    mask = (
        (ak.firsts(events.MuonGood.pt) > params["pt_leading_muon"])
        & OS
        & SF
        & (events.ll.mass > params["mll"]["low"])
        & (events.ll.mass < params["mll"]["high"])
        & (events.ll.deltaR > params["delta_r"])
    )

    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)


def PV_presel_cut(events, params, **kwargs):
    mask = abs(events.PV.z - events.GenVtx.z) < params["distance"]
    return ak.where(ak.is_none(mask), False, mask)


def get_custom_JetVetoMap_Mask(events, params, **kwargs):
    events_copy = copy.copy(events)
    events_copy["Jet"] = ak.with_field(
        events_copy["Jet"],
        events_copy[params["jet_type"]][params["pt_type"]],
        "pt",
    )

    mask = get_JetVetoMap_Mask(events_copy, params, **kwargs)

    return mask