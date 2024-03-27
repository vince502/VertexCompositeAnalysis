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
'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_30/src/jobConfigs/tmp/HIN-pPb816Summer16DR-00164_10.root',
#'file:output.root'
)
)

# =============== Other Statements =====================
# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
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
process.generalD0CandidatesNew.trkPtSumCut = cms.double(1.6)
process.generalD0CandidatesNew.trkEtaDiffCut = cms.double(2.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.6)
process.generalD0CandidatesNew.alphaCut = cms.double(2.0)
process.generalD0CandidatesNew.alpha2DCut = cms.double(2.0)
process.generalD0CandidatesNew.dPtCut = cms.double(1.9)

process.generalD0CandidatesNewWrongSign = process.generalD0CandidatesNew.clone(isWrongSign = cms.bool(True))

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDDCandidates_cff")
process.generalDDCandidatesNew = process.generalDDCandidates.clone()
process.generalDDCandidatesNew.trkPtSumCut = cms.double(0.0)
process.generalDDCandidatesNew.trkEtaDiffCut = cms.double(0.0)
process.generalDDCandidatesNew.tkNhitsCut = cms.int32(0)
process.generalDDCandidatesNew.tkPtErrCut = cms.double(0.0)
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

process.d0ana_mc.useAnyMVA = cms.bool(True)
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_mc.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")

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

process.generalDDCandidatesNew.d0Collection = cms.InputTag("d0selectorNewReduced:D0")

process.d0ana_mc_newreduced = process.d0ana_mc.clone()
process.d0ana_mc_newreduced.saveTree = True
process.d0ana_mc_newreduced.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNewReduced:D0")
process.d0ana_mc_newreduced.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")
process.d0ana_mc_newreduced.DCAValCollection = cms.InputTag("d0selectorMCNewReduced:DCAValuesNewD0")
process.d0ana_mc_newreduced.DCAErrCollection = cms.InputTag("d0selectorMCNewReduced:DCAErrorsNewD0")

process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorMCNewReduced * process.d0ana_mc_newreduced * process.generalDDCandidatesNew * process.ddana)

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

process.p = cms.Path(process.d0ana_seq2)

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

#process.outputPath = cms.EndPath(process.output)
#process.schedule.append(process.outputPath)
