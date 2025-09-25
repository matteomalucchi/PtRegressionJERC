import awkward as ak


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


def muon_selection_custom(events, params):
    lepton_flavour = "Muon"

    leptons = events[lepton_flavour]
    cuts = params.object_preselection[lepton_flavour]

    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]
    passes_dxy = abs(leptons.dxy) < cuts["dxy"]
    passes_dz = abs(leptons.dz) < cuts["dz"]
    passes_id = leptons[cuts["id"]] == True
    passes_iso = leptons.miniPFRelIso_all < cuts["iso"]
    passes_is_global = leptons.isGlobal == True
    passes_is_tracker = leptons.isTracker == True

    good_leptons = (
        passes_eta
        & passes_pt
        & passes_iso
        & passes_id
        & passes_dxy
        & passes_dz
        & passes_is_global
        & passes_is_tracker
    )

    return leptons[good_leptons]
