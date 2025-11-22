import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.kshortCandidates_cfi import KshortProducer

ChiCFromKshorts = cms.EDProducer(
    "ChiCResonanceProducer",
    resonanceCollection = cms.InputTag("KshortProducer", "Kshort"),
    applyMassWindow = cms.bool(True),
    requireUniqueTracks = cms.bool(True),
    # Vertex fitting configuration
    useVertexFitting = cms.bool(True),
    vertexRecoAlgorithm = cms.InputTag("offlinePrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    # Mass sigmas for kinematic fitting (for Ks: ~0.497 GeV with ~1% uncertainty)
    # These are used for the 4 tracks (2 pions from each Ks)
    resonanceMassSigmas = cms.vdouble(0.0013957018, 0.0013957018, 0.0013957018, 0.0013957018),  # Pion mass uncertainty (1%)
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
