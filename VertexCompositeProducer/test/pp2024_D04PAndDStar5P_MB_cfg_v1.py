import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.FastTimerService = cms.Service(
    "FastTimerService",
    printEventSummary=cms.untracked.bool(True),
    printRunSummary=cms.untracked.bool(True),
    printJobSummary=cms.untracked.bool(True),
    enableDQM=cms.untracked.bool(False),
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        '/store/data/Run2024J/PPRefZeroBiasPlusForward0/MINIAOD/PromptReco-v1/000/387/696/00000/0037fb37-713f-4df8-9668-a2ce4665a93c.root'
    ),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('141X_dataRun3_Express_v3')

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_PPRefZeroBias*',
]

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')
process.colEvtSel = cms.Sequence()

process.eventFilter_HM = cms.Sequence(process.hltFilter)
process.eventFilter_HM_step = cms.Path(process.eventFilter_HM)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

process.load('VertexCompositeAnalysis.VertexCompositeProducer.generalD04PCandidates_cff')
process.generalD04PCandidatesNew = process.generalD04PCandidates.clone()
process.generalD04PCandidatesNew.tkNhitsCut = cms.int32(6)
process.generalD04PCandidatesNew.tkPtCut = cms.double(0.6)
process.generalD04PCandidatesNew.tkPtErrCut = cms.double(0.08)
process.generalD04PCandidatesNew.tkEtaCut = cms.double(2.4)
process.generalD04PCandidatesNew.tkEtaDiffCut = cms.double(1.2)
process.generalD04PCandidatesNew.collinearityCut3D = cms.double(0.94)
process.generalD04PCandidatesNew.tkChi2Cut = cms.double(3.0)
process.generalD04PCandidatesNew.VtxChiProbCut = cms.double(0.01)
process.generalD04PCandidatesNew.vtxSignificance3DCut = cms.double(4.0)
process.generalD04PCandidatesNew.mPiKCutMin = cms.double(1.74)
process.generalD04PCandidatesNew.mPiKCutMax = cms.double(2.00)
process.generalD04PCandidatesNew.d0MassCut = cms.double(0.18)

process.load('VertexCompositeAnalysis.VertexCompositeProducer.generalDStar5PCandidates_cff')
process.generalDStar5PCandidatesNew = process.generalDStar5PCandidates.clone()
process.generalDStar5PCandidatesNew.d0Collection = cms.InputTag('generalD04PCandidatesNew:D04P')
process.generalDStar5PCandidatesNew.tkNhitsCut = cms.int32(3)
process.generalDStar5PCandidatesNew.tkPtCut = cms.double(0.45)
process.generalDStar5PCandidatesNew.tkChi2Cut = cms.double(3.0)
process.generalDStar5PCandidatesNew.dauLongImpactSigCut = cms.double(0.0)
process.generalDStar5PCandidatesNew.dauTransImpactSigCut = cms.double(0.0)
process.generalDStar5PCandidatesNew.dStarMassCut = cms.double(0.22)

process.d04prereco_step = cms.Path(process.eventFilter_HM * process.generalD04PCandidatesNew)

process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.d04panalyzer_tree_cff')
process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.dStar5panalyzer_tree_cff')
process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff')
process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.eventplaneanalyzer_cfi')

process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string('d04p_dstar5p_tree.root'),
)

process.d04pana_new = process.d04pana.clone()
process.d04pana_new.CompositeCollection = cms.untracked.InputTag('generalD04PCandidatesNew:D04P')
process.d04pana_new.MVACollection = cms.InputTag('generalD04PCandidatesNew:MVAValuesD04P')

process.dStar5pana_new = process.dStar5pana.clone()
process.dStar5pana_new.CompositeCollection = cms.untracked.InputTag('generalDStar5PCandidatesNew:DStar5P')
process.dStar5pana_new.MVACollection = cms.InputTag('generalDStar5PCandidatesNew:MVAValuesDStar5P')

process.dStar5PAna_step = cms.Path(
    process.eventFilter_HM
    * process.generalD04PCandidatesNew
    * process.generalDStar5PCandidatesNew
    * process.dStar5pana_new
)

process.d04pAna_step = cms.Path(
    process.eventFilter_HM
    * process.generalD04PCandidatesNew
    * process.d04pana_new
)

process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
process.eventinfoana.triggerPathNames = cms.untracked.vstring('HLT_PPRefZeroBias_v*')
process.eventinfoana.eventFilterNames = cms.untracked.vstring(
    'Flag_colEvtSel',
    'Flag_hfCoincFilter',
    'Flag_primaryVertexFilter',
)
process.eventinfoana.triggerFilterNames = cms.untracked.vstring()
process.eventinfoana.stageL1Trigger = cms.uint32(2)
process.pevt = cms.EndPath(process.eventinfoana)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d04pAna_step,
    process.dStar5PAna_step,
    process.pevt,
)

process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(
    process.eventFilter_HM * process.primaryVertexFilter * process.clusterCompatibilityFilter
)

eventFilterPaths = [process.Flag_colEvtSel, process.Flag_primaryVertexFilter]
for P in eventFilterPaths:
    process.schedule.insert(0, P)

changeToMiniAOD(process)
process.options.numberOfThreads = 1


process.MessageLogger.cerr.FwkReport.reportEvery = 1000