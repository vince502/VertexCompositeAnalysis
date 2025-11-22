import FWCore.ParameterSet.Config as cms

P4000Candidates = cms.EDProducer(
    "P4000Producer",
    jpsiCollection = cms.InputTag("DiMuonFromTracksProducer", "Jpsi"),
    phiCollection = cms.InputTag("DiKaonProducer", "DiKaon"),
    vertexRecoAlgorithm = cms.InputTag("offlinePrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    applyMassWindow = cms.bool(False),  # Store all candidates in range
    requireUniqueTracks = cms.bool(True),
    useVertexFitting = cms.bool(True),
    resonanceMassSigmas = cms.vdouble(0.003096916, 0.001019445),  # J/ψ and φ mass uncertainties
    states = cms.VPSet(
        cms.PSet(
            name = cms.string("P4000"),
            pdgId = cms.int32(9000443),  # Custom PDG ID for P(4000)
            mass = cms.double(4.0),
            massWindow = cms.double(0.5)  # ±500 MeV window
        )
    )
)
