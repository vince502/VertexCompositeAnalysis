import FWCore.ParameterSet.Config as cms

ChiCFlatNtuplizer = cms.EDAnalyzer(
    "ChiCFlatNtuplizer",
    treeName=cms.untracked.string("ChiCFlatNtuple"),
    primaryVertices=cms.InputTag("offlinePrimaryVertices"),
    sources=cms.VPSet(
        # To be filled in the job configuration.
    ),
)
