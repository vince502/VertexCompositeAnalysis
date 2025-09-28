import FWCore.ParameterSet.Config as cms

ChiCGenNtuplizer = cms.EDAnalyzer(
    'ChiCGenNtuplizer',
    genParticles=cms.InputTag('prunedGenParticles'),
    primaryVertices=cms.InputTag('offlinePrimaryVertices'),
    treeName=cms.untracked.string('ChiCGenNtuple')
)
