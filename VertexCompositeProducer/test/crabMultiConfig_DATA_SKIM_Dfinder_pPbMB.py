from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

from CRABClient.UserUtilities import config
config = config()

userName = "soohwan"
date = "2023Dec23"

config.section_("General")
config.General.workArea = 'crab_projects_knu_%s_v2' % (date)
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "pPbSkimAndTree2016_D0Both_BDT_ForLamC_cfg.py"

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
config.Data.runRange = '285479-285832'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
config.Site.whitelist = ['T2_US_Purdue', 'T2_US_Vanderbilt','T2_US_MIT']
#config.Site.blacklist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T3_KR_KNU'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print( "Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print( "Failed submitting task: %s" % (cle))

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {}
# MinBias 1 already ran in 06 March 2023
for i in range(1,21):
    dataMap[("PAMinimumBias"+str(i))] = { "PD": ("/PAMinimumBias"+str(i)+"/PARun2016C-PromptReco-v1/AOD"), "Units": 5, "Memory": 2500, "RunTime": 2400 }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = key+'_8TeVpPb2016_DmesonSkim_'+date+'_v4'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
#    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    #config.Data.outLFNDirBase = '/store/group/phys_heavyions/soohwan/Analysis/pPbMinimumBias/%s' % (config.General.requestName)
    config.Data.outLFNDirBase = '/store/user/soohwan/Analysis/%s/pPbMinimumBias/%s' % (date,config.General.requestName)

    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
