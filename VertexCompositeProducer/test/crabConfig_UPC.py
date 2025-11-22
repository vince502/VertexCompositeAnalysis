from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = "D04P_CMSSW_15_1_0_22Nov_2025_All_Forward2025_v1"
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = "Analysis"
#config.JobType.psetName = "UPC2025_ChiC_MB_cfg_v1.py"
config.JobType.psetName = "UPC2025_ChiC_Ultimate_cfg.py"
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2500         # request high memory machines.
#config.JobType.inputFiles=['CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374810.db']
#config.JobType.maxJobRuntimeMin = 2750    # request longer runtime, ~48 hours.

config.section_("Data")
#config.Data.inputDataset = '/HIForward4/HIRun2025A-PromptReco-v1/MINIAOD'
#config.Data.outputPrimaryDataset = 'ForwardCheck2025UPC'
config.Data.outputPrimaryDataset = 'ForwardCheck2024UPC'
config.Data.userInputFiles = open('forward.txt').readlines()
#config.Data.userInputFiles = open('forward2024.txt').readlines()
config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 16
config.Data.splitting = 'FileBased'
#config.Data.outLFNDirBase = '/store/user/soohwan/store/UPC2024/Charmonia/Chic/%s' % (config.General.requestName)
config.Data.outLFNDirBase = '/store/user/soohwan/store/UPC2025/Charmonia/Chic/%s' % (config.General.requestName)
config.Data.publication = False
config.Data.totalUnits = -1
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'

config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'
#config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_KR_KNU'
config.Site.whitelist = [ 'T2_US_*', 'T2_IT_*', 'T2_CH_*', 'T2_ES_*', 'T2_FR_*', 'T2_DE_*']
config.Site.ignoreGlobalBlacklist=True
