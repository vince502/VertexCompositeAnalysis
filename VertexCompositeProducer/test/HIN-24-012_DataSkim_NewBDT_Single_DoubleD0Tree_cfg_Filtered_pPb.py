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
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:PAMinimumBias19_AOD_PromptReco-v1_285975_000000_012C1E3-DCB5-E611-AE2F-02163E011ABE.root',
'/store/hidata/PARun2016C/PAMinimumBias8/AOD/PromptReco-v1/000/285/832/00001/0C1AABDC-B4B4-E611-8D01-02163E01460C.root',
)
)

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt').getVLuminosityBlockRange()

# =============== Other Statements =====================
# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '80X_dataRun2_v19'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # 'HLT_PAFullTracks_Multiplicity120_v*', # High multiplicity
    # 'HLT_PAFullTracks_Multiplicity150_v*', # High multiplicity
    'HLT_PAFullTracks_Multiplicity185_part*', # High multiplicity
    # 'HLT_PAFullTracks_Multiplicity250_v*', # High multiplicity
    'HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*', # Minimum bias
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
#process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)
#remove the default dz1p0 filter
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
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

process.generalD0CandidatesNew.tkPtSumCut = cms.double(1.6)
process.generalD0CandidatesNew.tkEtaDiffCut = cms.double(1.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(10)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.7)
process.generalD0CandidatesNew.alphaCut = cms.double(0.4)
process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)
process.generalD0CandidatesNew.dPtCut = cms.double(1.97)
process.generalD0CandidatesNew.d0MassCut = cms.double(0.20)
process.generalD0CandidatesNew.dAbsYCut = cms.double(1.01)
process.generalD0CandidatesNew.vtxSignificance3DCut = cms.double(1.00)
process.generalD0CandidatesNew.useAnyMVA = cms.bool(True)
process.generalD0CandidatesNew.mvaCut = cms.double(-1.0)
#process.generalD0CandidatesNew.GBRForestLabel = cms.string('D0InpPbXGB')
#process.generalD0CandidatesNew.GBRForestFileName = cms.string('GBRForestfile_XGBDT_PromptD0InpPb_default_MB_HardSoftQCD_wDauKine_v2_11Dec.root')
process.generalD0CandidatesNew.GBRForestLabel = cms.string('D0InpPb')
process.generalD0CandidatesNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')

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
    cms.string('ddana_data_tree.root')
    )

# set up selectors

process.ddana.PID = cms.untracked.int32(421)

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")

#process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v2.root')
#process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v2.root')
process.d0selectorMC = process.d0selectorMCBDTPreCut.clone()
process.d0selectorWS = process.d0selector.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
    MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.d0selectorMCNewReduced = process.d0selectorMC.clone()
#process.d0selectorMCNewReduced.GBRForestLabel = cms.string('D0InpPbXGB')
process.d0selectorMCNewReduced.GBRForestLabel = cms.string('D0InpPb')
#process.d0selectorMCNewReduced.GBRForestFileName = cms.string('GBRForestfile_XGBDT_PromptD0InpPb_default_MB_HardSoftQCD_wDauKine_v1.root')
process.d0selectorMCNewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.d0selectorMCNewReduced.DCAValCollection = cms.InputTag("generalD0CandidatesNew:DCAValuesD0")
process.d0selectorMCNewReduced.DCAErrCollection = cms.InputTag("generalD0CandidatesNew:DCAErrorsD0")
process.d0selectorMCNewReduced.trkPtMin = cms.untracked.double(0.7)
process.d0selectorMCNewReduced.trkPtSumMin = cms.untracked.double(0.0)
process.d0selectorMCNewReduced.trkEtaDiffMax = cms.untracked.double(1.0)
process.d0selectorMCNewReduced.trkNHitMin = cms.untracked.int32(10)
process.d0selectorMCNewReduced.cand3DPointingAngleMax = cms.untracked.double(1.0)
process.d0selectorMCNewReduced.cand2DPointingAngleMax = cms.untracked.double(1.0)
process.d0selectorMCNewReduced.candpTMin = cms.untracked.double(1.90)
process.d0selectorMCNewReduced.candYMin = cms.untracked.double(-1.11)
process.d0selectorMCNewReduced.candYMax = cms.untracked.double(1.11)
process.d0selectorMCNewReduced.mvaMin = cms.untracked.double(-10.0)


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

#process.generalDDCandidatesNew.d0Collection = cms.InputTag("d0selectorMCNewReduced:D0")
#process.generalDDCandidatesNew.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")

process.generalDDCandidatesNew.d0Collection = cms.InputTag("generalD0CandidatesNew:D0")
process.generalDDCandidatesNew.MVACollection = cms.InputTag("generalD0CandidatesNew:MVAValuesD0")
process.generalDDCandidatesNew.DCACollection = cms.InputTag("generalD0CandidatesNew:DCAValuesD0")
process.generalDDCandidatesNew.DCAErrCollection = cms.InputTag("generalD0CandidatesNew:DCAErrorsD0")
process.generalDDCandidatesNew.Angle2DCollection = cms.InputTag("generalD0CandidatesNew:Angle2D")
process.generalDDCandidatesNew.Angle3DCollection = cms.InputTag("generalD0CandidatesNew:Angle3D")


process.d0ana_newreduced = process.d0ana.clone()
process.d0ana_newreduced.saveTree = True
#process.d0ana_newreduced.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNewReduced:D0")
#process.d0ana_newreduced.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")
#process.d0ana_newreduced.DCAValCollection = cms.InputTag("d0selectorMCNewReduced:DCAValuesNewD0")
#process.d0ana_newreduced.DCAErrCollection = cms.InputTag("d0selectorMCNewReduced:DCAErrorsNewD0")

process.d0ana_newreduced.VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNew:D0")
process.d0ana_newreduced.MVACollection = cms.InputTag("generalD0CandidatesNew:MVAValuesD0")
process.d0ana_newreduced.DCAValCollection = cms.InputTag("generalD0CandidatesNew:DCAValuesD0")
process.d0ana_newreduced.DCAErrCollection = cms.InputTag("generalD0CandidatesNew:DCAErrorsD0")
process.d0ana_newreduced.Angle2D = cms.InputTag("generalD0CandidatesNew:Angle2D")
process.d0ana_newreduced.Angle3D = cms.InputTag("generalD0CandidatesNew:Angle3D")
process.d0ana_newreduced.MVACollection2= cms.InputTag("")
process.d0ana_newreduced.debug= False
process.d0ana_newreduced.onlyWantMatch= True

process.ddana_new = process.ddana.clone()
process.ddana_new.twoLayerDecay = cms.untracked.bool(True)
process.ddana_new.TrackCollection = cms.untracked.InputTag("generalTracks")
process.ddana_new.DCAValCollection = cms.InputTag("generalDDCandidatesNew:DCAValuesDD")
process.ddana_new.DCAErrCollection = cms.InputTag("generalDDCandidatesNew:DCAErrorsDD")
process.ddana_new.DCAValCollection1 = cms.InputTag("generalDDCandidatesNew:DCAValCollection1")
process.ddana_new.DCAErrCollection1 = cms.InputTag("generalDDCandidatesNew:DCAErrCollection1")
process.ddana_new.Angle2D1 = cms.InputTag("generalDDCandidatesNew:Angle2D1")
process.ddana_new.Angle3D1 = cms.InputTag("generalDDCandidatesNew:Angle3D1")
process.ddana_new.DCAValCollection2 = cms.InputTag("generalDDCandidatesNew:DCAValCollection2")
process.ddana_new.DCAErrCollection2 = cms.InputTag("generalDDCandidatesNew:DCAErrCollection2")
process.ddana_new.Angle2D2 = cms.InputTag("generalDDCandidatesNew:Angle2D2")
process.ddana_new.Angle3D2 = cms.InputTag("generalDDCandidatesNew:Angle3D2")
process.ddana_new.useAnyMVA = cms.bool(True)
# process.ddana_new.doGenMatching = cms.bool(True)
process.ddana_new.debug = cms.untracked.bool(False)
process.ddana_new.MVACollection = cms.InputTag("generalDDCandidatesNew:MVAValuesDD1")
process.ddana_new.MVACollection2= cms.InputTag("generalDDCandidatesNew:MVAValuesDD2")
process.ddana_new.onlyWantMatch = cms.untracked.bool(False)


#process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0ana_newreduced  )
process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0ana_newreduced * process.generalDDCandidatesNew * process.ddana_new )
#process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorMCNewReduced * process.d0ana_mc_newreduced  )
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

process.p = cms.Path(process.d0ana_seq2 )


# Add the Conversion tree

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.p,
    process.pevt,
)


# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter = cms.Path(process.eventFilter_HM * process.hfCoincFilter)
process.Flag_primaryVertexFilterPA = cms.Path(process.eventFilter_HM * process.primaryVertexFilterPA)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pPb8TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pPb8TeV_Gplus)
# follow the exactly same config of process.eventinfoana.eventFilterNames
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)

#process.output = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM"),
#    fileName = cms.untracked.string('output.root'),
#)
#
#process.outputPath = cms.EndPath(process.output)
#process.schedule.append(process.outputPath)

#process.options.numberOfThreads = cms.uint32(10)
