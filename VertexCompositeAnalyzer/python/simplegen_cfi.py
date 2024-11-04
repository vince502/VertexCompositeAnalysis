import FWCore.ParameterSet.Config as cms

genAna = cms.EDAnalyzer('GenParicleSimpleAnalyzer',
    GenParticleCollection = cms.untracked.InputTag("genParticles"),
    pdgIDs = cms.untracked.vint32([421, -421])
)