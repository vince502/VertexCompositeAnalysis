import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
#'/store/hidata/PARun2016C/PAMinimumBias19/AOD/PromptReco-v1/000/285/975/00000/0012C1E3-DCB5-E611-AE2F-02163E011ABE.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_10.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_11.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_12.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_13.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_2.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_3.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_4.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_5.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_6.root',
#'/store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/110000/682C6215-7F9A-E711-AB29-0025901D49AC.root',
#'/store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/110000/823C5ABC-839A-E711-A4B7-0CC47AC08BF8.root',
#'/store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/110000/D4B2A990-889A-E711-98C8-0CC47AC08BFA.root',
#'/store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/110000/88787D31-909A-E711-86AD-0025904C7B26.root',
#'/store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/110000/5C4D85CB-959A-E711-8E76-0CC47A7EEE80.root',
#'file:682C6215-7F9A-E711-AB29-0025901D49AC.root'
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/EvtGenMod/CMSSW_8_0_30/src/HIN-pPb816Summer16DR-00033.root'
#'file:HIN-pPb816Summer16DR-00033_68.root',
'/store/user/soohwan/Analysis/DDbarpPb/RECOAODSIM_MC_PromptD0_EvtGen_MOD_DDbar_PbPemb_25Aug2024_v2/PromptD0_Pythia8_SoftQCD_8TeVpPb_8_0_36_patch1/RECOAODSIM_MC_PromptD0_EvtGen_MOD_DDbar_PbPemb_25Aug2024_v2/240826_030023/0000/HIN-pPb816Summer16DR-00033_10.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_1.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_10.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_11.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_12.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_13.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_14.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_15.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_16.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_17.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_18.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_19.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_2.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_20.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_21.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_22.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_23.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_24.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_25.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_26.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_27.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_28.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_29.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_3.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_30.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_31.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_32.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_33.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_34.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_35.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_36.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_37.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_38.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_39.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_4.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_40.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_41.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_42.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_43.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_44.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_45.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_46.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_47.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_5.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_6.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_7.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_8.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_9.root',
)
)

# =============== Other Statements =====================
# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '80X_mcRun2_pA_v4'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping)


process.eventFilter_HM = cms.Sequence(
#    process.hltHM *
    process.colEvtSel
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
#process.generalD0CandidatesNew.trkPtSumCut = cms.double(1.6)
#process.generalD0CandidatesNew.trkEtaDiffCut = cms.double(2.0)
#process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
#process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
#process.generalD0CandidatesNew.tkPtCut = cms.double(0.6)
#process.generalD0CandidatesNew.alphaCut = cms.double(2.0)
#process.generalD0CandidatesNew.alpha2DCut = cms.double(2.0)
#process.generalD0CandidatesNew.dPtCut = cms.double(1.9)

process.generalD0CandidatesNew.trkPtSumCut = cms.double(0.0)
process.generalD0CandidatesNew.trkEtaDiffCut = cms.double(60.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(0)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(1.0)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.3)
process.generalD0CandidatesNew.alphaCut = cms.double(1000.0)
process.generalD0CandidatesNew.alpha2DCut = cms.double(1000.0)
process.generalD0CandidatesNew.dPtCut = cms.double(1.7)

process.generalD0CandidatesNewWrongSign = process.generalD0CandidatesNew.clone(isWrongSign = cms.bool(True))

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDDCandidates_cff")
process.generalDDCandidatesNew = process.generalDDCandidates.clone()
process.generalDDCandidatesNew.d0Collection= cms.InputTag("d0selectorMCNewReduced")
process.generalDDCandidatesNew.trkPtSumCut = cms.double(0.0)
process.generalDDCandidatesNew.trkEtaDiffCut = cms.double(999.0)
process.generalDDCandidatesNew.tkNhitsCut = cms.int32(0)
process.generalDDCandidatesNew.tkPtErrCut = cms.double(1.0)
process.generalDDCandidatesNew.tkPtCut = cms.double(0.0)
process.generalDDCandidatesNew.alphaCut = cms.double(999.0)
process.generalDDCandidatesNew.alpha2DCut = cms.double(999.0)
process.generalDDCandidatesNew.dPtCut = cms.double(0.0)


process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew)
process.d0rereco_wrongsign_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNewWrongSign )



# produce D0 trees
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cff")
#process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dStarselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.ddanalyzer_tree_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff")

process.TFileService = cms.Service("TFileService",
    fileName =
    cms.string('ddana_mc_tree.root')
    )

# set up selectors

process.ddana_mc.PID = cms.untracked.int32(421)

process.d0ana_mc.useAnyMVA = cms.bool(True)
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_mc.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")

#process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v2.root')
process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v2.root')
process.d0selectorMC = process.d0selectorMCBDTPreCut.clone()
process.d0selectorWS = process.d0selector.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
    MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.d0selectorMCNewReduced = process.d0selectorMC.clone()
process.d0selectorMCNewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.d0selectorMCNewReduced.DCAValCollection = cms.InputTag("generalD0CandidatesNew:DCAValuesD0")
process.d0selectorMCNewReduced.DCAErrCollection = cms.InputTag("generalD0CandidatesNew:DCAErrorsD0")
process.d0selectorMCNewReduced.trkPtMin = cms.untracked.double(0.3)
process.d0selectorMCNewReduced.trkPtSumMin = cms.untracked.double(0.0)
process.d0selectorMCNewReduced.trkEtaDiffMax = cms.untracked.double(100.0)
process.d0selectorMCNewReduced.trkNHitMin = cms.untracked.int32(0)
process.d0selectorMCNewReduced.cand3DPointingAngleMax = cms.untracked.double(100.0)
process.d0selectorMCNewReduced.cand2DPointingAngleMax = cms.untracked.double(100.0)
process.d0selectorMCNewReduced.candpTMin = cms.untracked.double(1.9)
process.d0selectorMCNewReduced.candYMin = cms.untracked.double(-1.1)
process.d0selectorMCNewReduced.candYMax = cms.untracked.double(1.1)
process.d0selectorMCNewReduced.mvaMin = cms.untracked.double(-1)


#process.d0selectorMCNewReduced.trkPtMin = cms.untracked.double(0.0)
#process.d0selectorMCNewReduced.trkPtSumMin = cms.untracked.double(0.0)
#process.d0selectorMCNewReduced.trkEtaDiffMax = cms.untracked.double(999.0)
#process.d0selectorMCNewReduced.trkNHitMin = cms.untracked.int32(0)
#process.d0selectorMCNewReduced.cand3DPointingAngleMax = cms.untracked.double(999.0)
#process.d0selectorMCNewReduced.cand2DPointingAngleMax = cms.untracked.double(999.0)
#process.d0selectorMCNewReduced.candpTMin = cms.untracked.double(0.0)
#process.d0selectorMCNewReduced.candYMin = cms.untracked.double(-5.0)
#process.d0selectorMCNewReduced.candYMax = cms.untracked.double(5.0)
#process.d0selectorMCNewReduced.mvaMin = cms.untracked.double(-1)

process.generalDDCandidatesNew.d0Collection = cms.InputTag("d0selectorMCNewReduced:D0")
process.generalDDCandidatesNew.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")


process.d0ana_mc_newreduced = process.d0ana_mc.clone()
process.d0ana_mc_newreduced.saveTree = True
process.d0ana_mc_newreduced.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNewReduced:D0")
process.d0ana_mc_newreduced.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")
process.d0ana_mc_newreduced.DCAValCollection = cms.InputTag("d0selectorMCNewReduced:DCAValuesNewD0")
process.d0ana_mc_newreduced.DCAErrCollection = cms.InputTag("d0selectorMCNewReduced:DCAErrorsNewD0")
process.d0ana_mc_newreduced.MVACollection2= cms.InputTag("")
process.d0ana_mc_newreduced.debug= False

process.ddana_new = process.ddana_mc.clone()
process.ddana_new.twoLayerDecay = cms.untracked.bool(True)
process.ddana_new.TrackCollection = cms.untracked.InputTag("generalTracks")
process.ddana_new.DCAValCollection = cms.InputTag("generalDDCandidatesNew:DCAValuesDD")
process.ddana_new.DCAErrCollection = cms.InputTag("generalDDCandidatesNew:DCAErrorsDD")
process.ddana_new.useAnyMVA = cms.bool(True)
# process.ddana_new.doGenMatching = cms.bool(True)
process.ddana_new.debug = cms.untracked.bool(False)
process.ddana_new.MVACollection = cms.InputTag("generalDDCandidatesNew:MVAValuesDD1")
process.ddana_new.MVACollection2= cms.InputTag("generalDDCandidatesNew:MVAValuesDD2")

process.ddana_old = process.ddana_mc_old.clone()
process.ddana_old.twoLayerDecay = cms.untracked.bool(True)
process.ddana_old.TrackCollection = cms.untracked.InputTag("generalTracks")
process.ddana_old.DCAValCollection = cms.InputTag("generalDDCandidatesNew:DCAValuesDD")
process.ddana_old.DCAErrCollection = cms.InputTag("generalDDCandidatesNew:DCAErrorsDD")
process.ddana_old.useAnyMVA = cms.bool(True)
# process.ddana_old.doGenMatching = cms.bool(True)
process.ddana_old.debug = cms.untracked.bool(False)
process.ddana_old.MVACollection = cms.InputTag("generalDDCandidatesNew:MVAValuesDD1")
process.ddana_old.MVACollection2= cms.InputTag("generalDDCandidatesNew:MVAValuesDD2")


process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")
process.d0ana_ntp_mc = process.d0ana_mc.clone()
process.d0ana_ntp_mc.useAnyMVA = cms.bool(True)
process.d0ana_ntp_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_ntp_mc.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")
process.d0ana_ntp_mc.saveTree = True
process.d0ana_ntp_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNewReduced:D0")
process.d0ana_ntp_mc.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")
process.d0ana_ntp_mc.DCAValCollection = cms.InputTag("d0selectorMCNewReduced:DCAValuesNewD0")
process.d0ana_ntp_mc.DCAErrCollection = cms.InputTag("d0selectorMCNewReduced:DCAErrorsNewD0")



process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorMCNewReduced * process.d0ana_mc_newreduced * process.d0ana_ntp_mc  )
#process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorMCNewReduced * process.d0ana_mc_newreduced) # * process.generalDDCandidatesNew * process.ddana_new )

# eventinfoana must be in EndPath, and process.eventinfoana.selectEvents must be the name of eventFilter_HM Path
process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
process.eventinfoana.triggerPathNames = cms.untracked.vstring(
    'HLT_PAFullTracks_Multiplicity120_v', # High multiplicity
    'HLT_PAFullTracks_Multiplicity150_v', # High multiplicity
    'HLT_PAFullTracks_Multiplicity185_part', # High multiplicity
    'HLT_PAFullTracks_Multiplicity250_v', # High multiplicity
    'HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part', # Minimum bias
    )
process.eventinfoana.triggerFilterNames = cms.untracked.vstring()
process.pevt = cms.EndPath(process.eventinfoana)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticles"),                                                                 
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(True),
                                   printStatus = cms.untracked.bool(False),
                                   printIndex = cms.untracked.bool(True),
                                   status = cms.untracked.vint32( 3 )
                                   )

process.p = cms.Path(process.d0ana_seq2 * process.printTree)


# Add the Conversion tree

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.p,
    process.pevt,
)

# # Add the event selection filters
# process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
# process.Flag_hfCoincFilter = cms.Path(process.eventFilter_HM * process.hfCoincFilter)
# process.Flag_primaryVertexFilterPA = cms.Path(process.eventFilter_HM * process.primaryVertexFilterPA)
# process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
# process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pPb8TeV_dz1p0)
# process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pPb8TeV_Gplus)
# # follow the exactly same config of process.eventinfoana.eventFilterNames
# eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
# for P in eventFilterPaths:
#     process.schedule.insert(0, P)

#process.output = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM"),
#    fileName = cms.untracked.string('output.root'),
#)
#
#process.outputPath = cms.EndPath(process.output)
#process.schedule.append(process.outputPath)

#process.options.numberOfThreads = cms.untracked.int32(8)
