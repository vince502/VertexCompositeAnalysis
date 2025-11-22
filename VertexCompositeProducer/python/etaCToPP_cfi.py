import FWCore.ParameterSet.Config as cms

EtaCToPP = cms.EDProducer(
    "ChiCTrackPairProducer",
    trackCollection = cms.InputTag("generalTracks"),
    daughterMass = cms.double(0.938272013),  # Proton mass
    minTrackPt = cms.double(0.5),
    maxTrackEta = cms.double(2.4),
    maxTrackNormalizedChi2 = cms.double(10.0),
    minTrackNHits = cms.int32(6),
    minTrackNPix = cms.int32(0),  # 0 = no cut (can be overridden)
    minPairPt = cms.double(1.0),
    applyMassWindow = cms.bool(True),
    requiredChargeProduct = cms.int32(-1),  # -1 for opposite charge (p+ p-)
    # Vertex fitting configuration
    useVertexFitting = cms.bool(True),
    vertexRecoAlgorithm = cms.InputTag("offlinePrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    daughterMassSigma = cms.double(0.00938272013),  # 1% uncertainty for proton mass
    states = cms.VPSet(
        cms.PSet(
            name = cms.string("EtaC"),
            pdgId = cms.int32(441),  # eta_c(1S) PDG ID
            mass = cms.double(2.9839),  # eta_c(1S) mass in GeV
            massWindow = cms.double(0.200)  # Wide window for initial studies
        )
    )
)
