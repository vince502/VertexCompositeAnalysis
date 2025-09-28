import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("CHICFLAT", eras.Run3_2023)

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = cms.string("141X_dataRun3_Express_v3")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(4),
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        "/store/data/Run2024J/PPRefHardProbes0/MINIAOD/PromptReco-v1/000/387/696/00000/0c419475-37b8-42d0-8fa1-9f576a68265d.root"
    ),
)

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string("chic_flat_ntuple.root"),
)

# --- Resonance producer: ChiC -> 4π ---
from VertexCompositeAnalysis.VertexCompositeProducer.chiCTo4Pi_cfi import ChiCTo4Pi as _ChiCTo4Pi

process.ChiCTo4Pi = _ChiCTo4Pi.clone(
    minTrackPt=cms.double(1.5),
    maxTrackEta=cms.double(2.4),
    maxTrackNormalizedChi2=cms.double(10.0),
    minTrackNHits=cms.int32(6),
    minCandidatePt=cms.double(3.5),
    minAcoplanarity=cms.double(0.6),
    maxSphericity=cms.double(0.35),
    maxCandidateAbsEta=cms.double(2.0),
    storeEventShape=cms.bool(True),
    applyMassWindow=cms.bool(True),
    states=cms.VPSet(
        cms.PSet(
            name=cms.string("ChiC0"),
            pdgId=cms.int32(10441),
            mass=cms.double(3.4147),
            massWindow=cms.double(0.025),
        ),
        cms.PSet(
            name=cms.string("ChiC2"),
            pdgId=cms.int32(445),
            mass=cms.double(3.5562),
            massWindow=cms.double(0.025),
        ),
    ),
)

# --- Flat ntuplizer ---
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCFlatNtuplizer_cfi import ChiCFlatNtuplizer as _ChiCFlatNtuplizer

process.ChiCFlatNtuplizer = _ChiCFlatNtuplizer.clone(
    treeName=cms.untracked.string("ChiCFlatNtuple"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    sources=cms.VPSet(
        cms.PSet(
            name=cms.string("ChiC0_4Pi"),
            pdgId=cms.int32(10441),
            collection=cms.InputTag("ChiCTo4Pi", "ChiC0"),
        ),
        cms.PSet(
            name=cms.string("ChiC2_4Pi"),
            pdgId=cms.int32(445),
            collection=cms.InputTag("ChiCTo4Pi", "ChiC2"),
        ),
    ),
)

# --- Analysis path ---
process.chicSequence = cms.Sequence(process.ChiCTo4Pi)
process.analysisPath = cms.Path(process.chicSequence * process.ChiCFlatNtuplizer)

process.schedule = cms.Schedule(process.analysisPath)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

changeToMiniAOD(process)
