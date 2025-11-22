import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('ANASKIM', eras.Run3_2025)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(1),
)
process.FastTimerService = cms.Service(
    'FastTimerService',
    printEventSummary=cms.untracked.bool(True),
    printRunSummary=cms.untracked.bool(True),
    printJobSummary=cms.untracked.bool(True),
    enableDQM=cms.untracked.bool(False),
)

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
#'/store/hidata/OORun2025/IonPhysics0/MINIAOD/PromptReco-v1/000/394/075/00000/09db905b-c8ac-4e9e-9d6d-2be7f844a12b.root'
#'file:04e18742-3308-45a5-b0d6-560741bec33f.root',
        # '/store/data/Run2024J/PPRefZeroBiasPlusForward0/MINIAOD/PromptReco-v1/000/387/696/00000/0037fb37-713f-4df8-9668-a2ce4665a93c.root'
'/store/hidata/HIRun2025A/HIForward0/MINIAOD/PromptReco-v1/000/399/540/00000/491ce449-be44-4fe8-a337-f85c90499ea9.root',
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
]

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')
process.colEvtSel = cms.Sequence()

# Track count filter - configurable cut on number of selected tracks
from VertexCompositeAnalysis.VertexCompositeProducer.trackCountFilter_cfi import trackCountFilter
process.trackCountFilter = trackCountFilter.clone()
# Configure track selection cuts (matching ChiC producers)
process.trackCountFilter.trackCollection = cms.InputTag('generalTracks')
process.trackCountFilter.minTrackPt = cms.double(0.1)  # Match ChiCTo4Pi/ChiCTo2Ka
process.trackCountFilter.maxTrackEta = cms.double(2.4)  # Match ChiC2To4K/ChiCTo2Ka
process.trackCountFilter.maxTrackNormalizedChi2 = cms.double(1000.0)
process.trackCountFilter.minTrackNHits = cms.int32(0)  # Match ChiCTo4Pi/ChiC2To4K/ChiCTo2Ka
# Configure track count cuts
process.trackCountFilter.minNSelectedTracks = cms.int32(4)  # Minimum tracks required (>=)
process.trackCountFilter.maxNSelectedTracks = cms.int32(6)  # Maximum tracks (<=)

process.eventFilter_HM = cms.Sequence(process.hltFilter * process.trackCountFilter)
process.eventFilter_HM_step = cms.Path(process.eventFilter_HM)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# --- D0 and B producers ---
process.load('VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff')
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(0)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.3)
process.generalD0CandidatesNew.tkEtaCut = cms.double(2.4)
process.generalD0CandidatesNew.tkChi2Cut = cms.double(10.0)
process.generalD0CandidatesNew.d0MassCut = cms.double(0.15)

process.load('VertexCompositeAnalysis.VertexCompositeProducer.generalBCandidates_cff')
process.generalBCandidatesNew = process.generalBCandidates.clone()
process.generalBCandidatesNew.d0RecoAlgorithm = cms.InputTag('generalD0CandidatesNew', 'D0')
process.generalBCandidatesNew.batTkPtCut = cms.double(0.3)
process.generalBCandidatesNew.batTkEtaCut = cms.double(2.4)
process.generalBCandidatesNew.batTkChi2Cut = cms.double(10.0)
process.generalBCandidatesNew.batTkNhitsCut = cms.int32(0)

# --- ChiC producers (4pi only) ---
from VertexCompositeAnalysis.VertexCompositeProducer.chiCTo4Pi_cfi import ChiCTo4Pi as _ChiCTo4Pi
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCFlatNtuplizer_cfi import ChiCFlatNtuplizer as _ChiCFlatNtuplizer

process.ChiCTo4Pi = _ChiCTo4Pi.clone()
# Expose key ChiCTo4Pi selections
process.ChiCTo4Pi.minTrackPt = cms.double(0.1)
process.ChiCTo4Pi.maxTrackEta = cms.double(2.4)
process.ChiCTo4Pi.maxTrackNormalizedChi2 = cms.double(999.0)
process.ChiCTo4Pi.minTrackNHits = cms.int32(0)
process.ChiCTo4Pi.minCandidatePt = cms.double(0)
process.ChiCTo4Pi.minAcoplanarity = cms.double(0.0)
process.ChiCTo4Pi.maxSphericity = cms.double(1.0)
process.ChiCTo4Pi.maxCandidateAbsEta = cms.double(100.0)
process.ChiCTo4Pi.storeEventShape = cms.bool(True)
process.ChiCTo4Pi.applyMassWindow = cms.bool(True)
# Enable vertex fitting to get quality metrics
process.ChiCTo4Pi.useVertexFitting = cms.bool(True)
process.ChiCTo4Pi.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.ChiCTo4Pi.beamSpot = cms.InputTag('offlineBeamSpot')
# Single state since we don't care about mass cuts - all 4-pion candidates go to one collection
# If you need separate ChiC0/ChiC2 collections or different mass windows, use multiple states
process.ChiCTo4Pi.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC'),
        pdgId=cms.int32(445),  # ChiC2 PDG ID (or use 10441 for ChiC0, doesn't matter if not using mass cuts)
        mass=cms.double(3.5),  # Not used if applyMassWindow=False
        massWindow=cms.double(10.0)  # Wide window, not used if applyMassWindow=False
    )
)

# Unified ChiC flat ntuple writer (4pi only) - comprehensive information for offline analysis
# Single collection since we're using single state
process.ChiCFlatNtuplizer = _ChiCFlatNtuplizer.clone(
    treeName=cms.untracked.string('ChiCFlatNtuple'),
    primaryVertices=cms.InputTag('offlinePrimaryVertices'),
    sources=cms.VPSet(
        cms.PSet(
            name=cms.string('ChiC_4Pi'),
            pdgId=cms.int32(445),  # Match the pdgId in states
            collection=cms.InputTag('ChiCTo4Pi', 'ChiC')
        )
    )
)

# --- Analysis paths ---
# D0 and B reconstruction
process.d0_step = cms.Path(
    process.eventFilter_HM * process.generalD0CandidatesNew
)

process.b_step = cms.Path(
    process.eventFilter_HM * process.generalD0CandidatesNew * process.generalBCandidatesNew
)

# ChiC to 4pi only
process.chic4Pi_step = cms.Path(
    process.eventFilter_HM * process.ChiCTo4Pi
)

process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff')
process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string('chic_combination_tree.root'),
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
    #process.d0_step,
    #process.b_step,
    process.chic4Pi_step,
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
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.outCustom = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("myOutput.root"),
    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM")  # Keep everything
)

# Output path
process.chicNtuple = cms.EndPath(process.ChiCFlatNtuplizer)
#process.outpathcustom = cms.EndPath(process.outCustom)
#process.schedule.append(process.outpathcustom)
process.schedule.append(process.chicNtuple)
