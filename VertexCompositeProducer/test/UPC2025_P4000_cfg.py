import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('ANASKIM', eras.Run3_2025_UPC)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(1),
)

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        '/store/hidata/HIRun2025A/HIForward0/MINIAOD/PromptReco-v1/000/399/655/00000/04672a9e-7153-4e94-8603-98e3eb08fa58.root'
    ),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('151X_dataRun3_Prompt_v1')

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
    '*',
]

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')
process.colEvtSel = cms.Sequence()

# Track count filter
from VertexCompositeAnalysis.VertexCompositeProducer.trackCountFilter_cfi import trackCountFilter
process.trackCountFilter = trackCountFilter.clone()
process.trackCountFilter.trackCollection = cms.InputTag('generalTracks')
process.trackCountFilter.minTrackPt = cms.double(0.1)
process.trackCountFilter.maxTrackEta = cms.double(2.4)
process.trackCountFilter.maxTrackNormalizedChi2 = cms.double(10.0)
process.trackCountFilter.minTrackNHits = cms.int32(5)
process.trackCountFilter.minTrackNPix = cms.int32(3)
process.trackCountFilter.minNSelectedTracks = cms.int32(2)
process.trackCountFilter.maxNSelectedTracks = cms.int32(10)

process.eventFilter_HM = cms.Sequence(process.hltFilter * process.trackCountFilter)
process.eventFilter_HM_step = cms.Path(process.eventFilter_HM)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# ============================================================================
# P(4000) -> J/ψ(μ+μ-) + φ(K+K-) Analysis Chain
# ============================================================================

# 1. Di-Muon Producer (J/ψ from tracks with muon hypothesis)
from VertexCompositeAnalysis.VertexCompositeProducer.diMuonFromTracksCandidates_cfi import DiMuonFromTracksCandidates
process.DiMuonFromTracksProducer = DiMuonFromTracksCandidates.clone()
process.DiMuonFromTracksProducer.trackRecoAlgorithm = cms.InputTag('generalTracks')  # Will be replaced by changeToMiniAOD
process.DiMuonFromTracksProducer.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')  # Will be replaced by changeToMiniAOD

# 2. Di-Kaon Producer (φ from tracks with kaon hypothesis)
from VertexCompositeAnalysis.VertexCompositeProducer.diKaonCandidates_cfi import DiKaonProducer as _DiKaonProducer
process.DiKaonProducer = _DiKaonProducer.clone()
process.DiKaonProducer.trackRecoAlgorithm = cms.InputTag('generalTracks')  # Will be replaced by changeToMiniAOD
process.DiKaonProducer.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')  # Will be replaced by changeToMiniAOD
process.DiKaonProducer.tkChi2Cut = cms.double(10.0)
process.DiKaonProducer.tkNhitsCut = cms.int32(5)
process.DiKaonProducer.tkPtCut = cms.double(0.1)
process.DiKaonProducer.tkEtaCut = cms.double(2.4)
process.DiKaonProducer.phiMassCut = cms.double(1.0)  # Loose cut
process.DiKaonProducer.mKKCutMin = cms.double(0.2)
process.DiKaonProducer.mKKCutMax = cms.double(2.0)

# 3. P(4000) Producer (combines J/ψ and φ)
from VertexCompositeAnalysis.VertexCompositeProducer.p4000Candidates_cfi import P4000Candidates
process.P4000Candidates = P4000Candidates.clone()
process.P4000Candidates.jpsiCollection = cms.InputTag('DiMuonFromTracksProducer', 'Jpsi')
process.P4000Candidates.phiCollection = cms.InputTag('DiKaonProducer', 'DiKaon')
process.P4000Candidates.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')  # Will be replaced by changeToMiniAOD
process.P4000Candidates.beamSpot = cms.InputTag('offlineBeamSpot')
process.P4000Candidates.applyMassWindow = cms.bool(False)
process.P4000Candidates.requireUniqueTracks = cms.bool(True)
process.P4000Candidates.useVertexFitting = cms.bool(True)

# 4. Flat Ntuplizer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.p4000FlatNtuplizer_cfi import P4000FlatNtuplizer
process.P4000FlatNtuplizer = P4000FlatNtuplizer.clone()
process.P4000FlatNtuplizer.p4000Collection = cms.InputTag('P4000Candidates', 'P4000')
process.P4000FlatNtuplizer.primaryVertices = cms.InputTag('offlinePrimaryVertices')  # Will be replaced by changeToMiniAOD

# ============================================================================
# Analysis Paths
# ============================================================================
process.p4000_step = cms.Path(
    process.eventFilter_HM *
    process.DiMuonFromTracksProducer *
    process.DiKaonProducer *
    process.P4000Candidates
)

process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff')
process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string('p4000_tree.root'),
)

process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
process.eventinfoana.triggerPathNames = cms.untracked.vstring('HLT_*')
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
    process.p4000_step,
    process.pevt,
)

process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(
    #process.eventFilter_HM * process.primaryVertexFilter * process.clusterCompatibilityFilter
    process.eventFilter_HM * process.primaryVertexFilter
)

eventFilterPaths = [process.Flag_colEvtSel, process.Flag_primaryVertexFilter]
for P in eventFilterPaths:
    process.schedule.insert(0, P)

changeToMiniAOD(process)
# After changeToMiniAOD, ensure ntuplizer uses the unpacked vertices
process.P4000FlatNtuplizer.primaryVertices = cms.InputTag('unpackedTracksAndVertices')
process.options.numberOfThreads = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Output path for flat ntuple
process.p4000Ntuple = cms.EndPath(process.P4000FlatNtuplizer)
process.schedule.append(process.p4000Ntuple)
