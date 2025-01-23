import FWCore.ParameterSet.Config as cms

PVMultiFSana = cms.EDAnalyzer('PrimaryVertexAnalyzerMultiFS',
    VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
    TrackCollection = cms.untracked.InputTag("generalTracks"),
    VertexCompositeCollection = cms.untracked.InputTag("VertexCompositeCollection"),
    doGen = cms.bool(False),
    doReco = cms.bool(False),
)