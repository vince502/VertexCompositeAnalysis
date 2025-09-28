import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.diKaonCandidates_cfi import DiKaonProducer

ChiCFromDiKaons = cms.EDProducer(
    "ChiCResonanceProducer",
    resonanceCollection = cms.InputTag("DiKaonProducer", "DiKaon"),
    applyMassWindow = cms.bool(True),
    requireUniqueTracks = cms.bool(True),
    states = cms.VPSet(
        cms.PSet(
            name = cms.string("ChiC0"),
            pdgId = cms.int32(10441),
            mass = cms.double(3.4147),
            massWindow = cms.double(0.060)
        ),
        cms.PSet(
            name = cms.string("ChiC2"),
            pdgId = cms.int32(445),
            mass = cms.double(3.5562),
            massWindow = cms.double(0.060)
        )
    )
)

ChiCFromDiKaonsSequence = cms.Sequence(
    DiKaonProducer * ChiCFromDiKaons
)
