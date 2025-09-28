import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.kshortCandidates_cfi import KshortProducer

ChiCFromKshorts = cms.EDProducer(
    "ChiCResonanceProducer",
    resonanceCollection = cms.InputTag("KshortProducer", "Kshort"),
    applyMassWindow = cms.bool(True),
    requireUniqueTracks = cms.bool(True),
    states = cms.VPSet(
        cms.PSet(
            name = cms.string("ChiC0"),
            pdgId = cms.int32(10441),
            mass = cms.double(3.4147),
            massWindow = cms.double(0.050)
        ),
        cms.PSet(
            name = cms.string("ChiC2"),
            pdgId = cms.int32(445),
            mass = cms.double(3.5562),
            massWindow = cms.double(0.050)
        )
    )
)

ChiCFromKshortsSequence = cms.Sequence(
    KshortProducer * ChiCFromKshorts
)
