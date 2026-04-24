from pocket_coffea.parameters.histograms import HistConf, Axis
import numpy as np

from mc_truth_ptreg_jerc.params.binning import *

# Define rho bins (density of pileup, typically ranges 0-100)
rho_bins = np.linspace(0, 100, 100)

def get_variables_dict(cuts_names_eta, cuts_names_reco_eta, cuts_names_eta_neutrino, mc_truth=True):
    """
    Generate histograms for JER analysis with N-dimensional binning.
    
    Histograms are organized as:
    - rho_pu: (PU, rho) - for rho vs pileup linear fit
    - Response histograms: (PU, eta, pt, response_variable) - for resolution computation
    - Each response type has its own histogram
    
    Returns dictionaries keyed like 'h_rho_pu', 'h_ResponseJEC', etc. when used.
    """
    variables_dict = {}
    
    # ===== RHO VS PILEUP (for linear fit) =====
    # This histogram has Mean storage to track average rho in pileup bins
    variables_dict["rho_pu"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="Rho",
                field="fixedGridRhoFastjetAll",
                bins=rho_bins,
                label="Rho_fixedGridRhoFastjetAll",
            ),
        ]
    )
    
    # ===== STANDARD JET PT VARIABLES (MatchedJets) =====
    # (PU, eta, pt, reco_pt) - 4D histogram
    # JetPtJEC histograms 
    
    # variables_dict["JetPtJEC"] = HistConf(
    #     [
            
    
    
    
    # ===== STANDARD RESPONSE VARIABLES (MatchedJets) =====
    # (PU, eta, pt, response) - 4D histogram
    
    # ResponseJEC
    variables_dict["ResponseJEC"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJets",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJets_eta",
            ),
            Axis(
                coll="MatchedJets",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJets_pt",
            ),
            Axis(
                coll="MatchedJets",
                field="ResponseJEC",
                bins=response_bins,
                label="MatchedJets_ResponseJEC",
            ),
        ]
    )
    
    # ResponseRaw
    variables_dict["ResponseRaw"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJets",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJets_eta",
            ),
            Axis(
                coll="MatchedJets",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJets_pt",
            ),
            Axis(
                coll="MatchedJets",
                field="ResponseRaw",
                bins=response_bins,
                label="MatchedJets_ResponseRaw",
            ),
        ]
    )
    
    # ResponsePNetReg
    variables_dict["ResponsePNetReg"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJets",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJets_eta",
            ),
            Axis(
                coll="MatchedJets",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJets_pt",
            ),
            Axis(
                coll="MatchedJets",
                field="ResponsePNetReg",
                bins=response_bins,
                label="MatchedJets_ResponsePNetReg",
            ),
        ]
    )
    
    # ResponseUparTReg
    variables_dict["ResponseUparTReg"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJets",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJets_eta",
            ),
            Axis(
                coll="MatchedJets",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJets_pt",
            ),
            Axis(
                coll="MatchedJets",
                field="ResponseUparTReg",
                bins=response_bins,
                label="MatchedJets_ResponseUparTReg",
            ),
        ]
    )
    
    
    # ===== NEUTRINO RESPONSE VARIABLES (MatchedJetsNeutrino) =====
    # (PU, eta, pt, response) - 4D histogram
    
    # ResponseJEC (Neutrino)
    variables_dict["ResponseJEC_neutrino"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJetsNeutrino_eta",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJetsNeutrino_pt",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="ResponseJEC",
                bins=response_bins,
                label="MatchedJetsNeutrino_ResponseJEC",
            ),
        ]
    )
    
    # ResponseRaw (Neutrino)
    variables_dict["ResponseRaw_neutrino"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJetsNeutrino_eta",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJetsNeutrino_pt",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="ResponseRaw",
                bins=response_bins,
                label="MatchedJetsNeutrino_ResponseRaw",
            ),
        ]
    )
    
    # ResponsePNetRegNeutrino
    variables_dict["ResponsePNetRegNeutrino"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJetsNeutrino_eta",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJetsNeutrino_pt",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="ResponsePNetRegNeutrino",
                bins=response_bins,
                label="MatchedJetsNeutrino_ResponsePNetRegNeutrino",
            ),
        ]
    )
    
    # ResponseUparTRegNeutrino
    variables_dict["ResponseUparTRegNeutrino"] = HistConf(
        [
            Axis(
                coll="Pileup",
                field="nPU",
                bins=PU_bins,
                label="Pileup_nPU",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="eta",
                bins=jet_eta_bins,
                label="MatchedJetsNeutrino_eta",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="pt",
                bins=jet_pt_bins,
                label="MatchedJetsNeutrino_pt",
            ),
            Axis(
                coll="MatchedJetsNeutrino",
                field="ResponseUparTRegNeutrino",
                bins=response_bins,
                label="MatchedJetsNeutrino_ResponseUparTRegNeutrino",
            ),
        ]
    )
    
    return variables_dict