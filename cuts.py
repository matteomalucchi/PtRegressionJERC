from pocket_coffea.lib.cut_definition import Cut
import awkward as ak
import os
import custom_cut_functions as cf


def get_ptbin(pt_low, pt_high, name=None):
    if name == None:
        name = f"pt{pt_low}to{pt_high}"
    return Cut(
        name=name,
        params={"pt_low": pt_low, "pt_high": pt_high},
        function=cf.ptbin,
        collection="MatchedJets",
    )


def get_etabin(eta_low, eta_high, name=None):
    if name == None:
        name = f"MatchedJets_eta{eta_low}to{eta_high}"
    return Cut(
        name=name,
        params={"eta_low": eta_low, "eta_high": eta_high},
        function=cf.etabin,
        collection=("MatchedJets"),
    )


def get_etabin_neutrino(eta_low, eta_high, name=None):
    if name == None:
        name = f"MatchedJetsNeutrino_eta{eta_low}to{eta_high}"
    return Cut(
        name=name,
        params={"eta_low": eta_low, "eta_high": eta_high},
        function=cf.etabin_neutrino,
        collection="MatchedJetsNeutrino",
    )


def get_reco_etabin(eta_low, eta_high, name=None):
    if name == None:
        name = f"MatchedJets_reco_eta{eta_low}to{eta_high}"
    return Cut(
        name=name,
        params={"eta_low": eta_low, "eta_high": eta_high},
        function=cf.reco_etabin,
        collection=("MatchedJets"),
    )


def get_reco_neutrino_etabin(eta_low, eta_high, name=None):
    if name == None:
        name = f"MatchedJetsNeutrino_reco_eta{eta_low}to{eta_high}"
    return Cut(
        name=name,
        params={"eta_low": eta_low, "eta_high": eta_high},
        function=cf.reco_neutrino_etabin,
        collection=("MatchedJetsNeutrino"),
    )


def get_reco_neutrino_abs_etabin(eta_low, eta_high, name=None):
    if name == None:
        name = f"MatchedJetsNeutrino_reco_abs_eta{eta_low}to{eta_high}"
    return Cut(
        name=name,
        params={"eta_low": eta_low, "eta_high": eta_high},
        function=cf.reco_neutrino_abs_etabin,
        collection=("MatchedJetsNeutrino"),
    )

PV_presel = Cut(
    name="PV_presel",
    params={
        "distance": 0.2
    },
    function=cf.PV_presel_cuts,
)