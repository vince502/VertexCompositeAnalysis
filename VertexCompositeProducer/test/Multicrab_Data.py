from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
# from httplib import HTTPException
from http.client import HTTPException

from CRABClient.UserUtilities import config
config = config()

userName = "junseok"
date = "20241010"

config.section_("General")
config.General.workArea = 'crab_projects/'+date
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374810.db']
config.JobType.psetName = 'PbPb2023_D0BothAndDStar_MB_cfg_v1.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
# config.JobType.pyCfgParams = [ 'numThreads=1', 'resonance=Z', 'isFullAOD=False', 'isMC=False', 'globalTag=132X_dataRun3_Prompt_v4', 'era=Run2023HI', 'includeJets=False' ]
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
#config.Data.runRange = '326381-327564'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_US_Purdue', 'T2_US_MIT']
#config.Site.blacklist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print ("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print ("Failed submitting task: %s" % (cle))

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {}

for i in range(10,32):
    dataMap[("HIPhysicsRawPrime"+str(i))] = { "PD": ("/HIPhysicsRawPrime"+str(i)+"/HIRun2023A-PromptReco-v2/MINIAOD"), "Units": 20, "Memory": 2500, "RunTime": 800 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'DStarAnalysis_PbPb2023_DstarToKpipi_CMSSW_13_2_11_'+key+'_HIRun2023A-PromptReco'+date
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    # config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/Run2023HI/MINIAOD/%s/%s' % (userName, date, config.General.requestName)
    print("Submitting CRAB job for: "+val["PD"])
    submit(config)