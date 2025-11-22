import FWCore.ParameterSet.Config as cms

P4000FlatNtuplizer = cms.EDAnalyzer(
    "P4000FlatNtuplizer",
    p4000Collection = cms.InputTag("P4000Candidates", "P4000"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices")
)
