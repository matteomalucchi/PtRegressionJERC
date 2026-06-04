import os


localdir = os.path.dirname(os.path.abspath(__file__))


corr_files_pnetreg = {
    "2016_PreVFP": "Summer20UL16APV_V1_MC_L2Relative_AK4PFPNet.txt",
    "2016_PostVFP": "Summer20UL16_V1_MC_L2Relative_AK4PFPNet.txt",
    "2017": "Summer20UL17_V1_MC_L2Relative_AK4PFPNet.txt",
    "2017": "Summer20UL17_V1_MC_L2Relative_AK4PFPNet.txt",
    "2018": "Summer20UL18_V1_MC_L2Relative_AK4PFPNet.txt",
    "2022_preEE": "Summer22Run3_V3_MC_L2Relative_AK4PFPNet.txt",
    "2022_postEE": "Summer22EERun3_V3_MC_L2Relative_AK4PFPNet.txt",
    "2023_preBPix": "Summer23Run3_V3_MC_L2Relative_AK4PFPNet.txt",
    "2023_postBPix": "Summer23BPixRun3_V3_MC_L2Relative_AK4PFPNet.txt",
    "2024": "Summer24Run3_V3_MC_L2Relative_AK4PFPNet.txt",
    "2025": "Winter25Run3_V3_MC_L2Relative_AK4PFPNet.txt",
}

corr_files_upartreg = {
    "2016_PreVFP": "Summer20UL16APV_V1_MC_L2Relative_AK4PFUparT.txt",
    "2016_PostVFP": "Summer20UL16_V1_MC_L2Relative_AK4PFUparT.txt",
    "2017": "Summer20UL17_V1_MC_L2Relative_AK4PFUparT.txt",
    "2018": "Summer20UL18_V1_MC_L2Relative_AK4PFUparT.txt",
    "2022_preEE": "Summer22Run3_V3_MC_L2Relative_AK4PFUparT.txt",
    "2022_postEE": "Summer22EERun3_V3_MC_L2Relative_AK4PFUparT.txt",
    "2023_preBPix": "Summer23Run3_V3_MC_L2Relative_AK4PFUparT.txt",
    "2023_postBPix": "Summer23BPixRun3_V3_MC_L2Relative_AK4PFUparT.txt",
    "2024": "Summer24Run3_V3_MC_L2Relative_AK4PFUparT.txt",
    "2025": "Winter25Run3_V3_MC_L2Relative_AK4PFUparT.txt",
}

corr_files_pnetreg_neutrino = {
    "2016_PreVFP": "Summer20UL16APV_V1_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2016_PostVFP": "Summer20UL16_V1_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2017": "Summer20UL17_V1_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2018": "Summer20UL18_V1_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2022_preEE": "Summer22Run3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2022_postEE": "Summer22EERun3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2023_preBPix": "Summer23Run3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2023_postBPix": "Summer23BPixRun3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2024": "Summer24Run3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
    "2025": "Winter25Run3_V3_MC_L2Relative_AK4PFPNetPlusNeutrino.txt",
}

corr_files_upartreg_neutrino = {
    "2016_PreVFP": "Summer20UL16APV_V1_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2016_PostVFP": "Summer20UL16_V1_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2017": "Summer20UL17_V1_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2018": "Summer20UL18_V1_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2022_preEE": "Summer22Run3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2022_postEE": "Summer22EERun3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2023_preBPix": "Summer23Run3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2023_postBPix": "Summer23BPixRun3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2024": "Summer24Run3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
    "2025": "Winter25Run3_V3_MC_L2Relative_AK4PFUparTPlusNeutrino.txt",
}

corr_files = {
    "2016_PreVFP": "Summer20UL16APV_V1_MC_L2Relative_AK4PFPuppi.txt",
    "2016_PostVFP": "Summer20UL16_V1_MC_L2Relative_AK4PFPuppi.txt",
    "2017": "Summer20UL17_V1_MC_L2Relative_AK4PFPuppi.txt",
    "2018": "Summer20UL18_V1_MC_L2Relative_AK4PFPuppi.txt",
    "2022_preEE": "Summer22Run3_V1_MC_L2Relative_AK4PUPPI.txt",
    "2022_postEE": "Summer22EEVetoRun3_V1_MC_L2Relative_AK4PUPPI.txt",
    "2023_preBPix": "Summer23Run3_V1_MC_L2Relative_AK4PUPPI.txt",
    "2023_postBPix": "Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI.txt",
    "2024": "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi.txt",
    "2025": "Winter25Prompt25_V1_MC_L2Relative_AK4PFPuppi.txt",
}

def get_corr_file(base_path):
    
    corr_files_pnetreg_full_path = {}
    corr_files_upartreg_full_path = {}
    corr_files_pnetreg_neutrino_full_path = {}
    corr_files_upartreg_neutrino_full_path = {}
    corr_files_full_path = {}
    
    years = list(corr_files.keys())

    for year in years:
        corr_files_pnetreg_full_path[year] = os.path.join(base_path, corr_files_pnetreg[year])
        corr_files_upartreg_full_path[year] = os.path.join(base_path, corr_files_upartreg[year])
        corr_files_pnetreg_neutrino_full_path[year] = os.path.join(base_path, corr_files_pnetreg_neutrino[year])
        corr_files_upartreg_neutrino_full_path[year] = os.path.join(base_path, corr_files_upartreg_neutrino[year])
        corr_files_full_path[year] = os.path.join(base_path, corr_files[year])
        
    return corr_files_pnetreg_full_path, corr_files_upartreg_full_path, corr_files_pnetreg_neutrino_full_path, corr_files_upartreg_neutrino_full_path, corr_files_full_path