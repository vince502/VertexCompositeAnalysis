import FWCore.ParameterSet.Config as cms

ChiCTo2Ka = cms.EDProducer("ChiCTrackPairProducer")

# Track selection cuts
ChiCTo2Ka.trackCollection = cms.InputTag("generalTracks")
ChiCTo2Ka.daughterMass = cms.double(0.493677)
ChiCTo2Ka.minTrackPt = cms.double(0.6)
ChiCTo2Ka.maxTrackEta = cms.double(2.4)
ChiCTo2Ka.maxTrackNormalizedChi2 = cms.double(10.0)
ChiCTo2Ka.minTrackNHits = cms.int32(6)

# Pair-level cuts
ChiCTo2Ka.minPairPt = cms.double(1.2)
ChiCTo2Ka.applyMassWindow = cms.bool(True)
ChiCTo2Ka.requiredChargeProduct = cms.int32(-1)  # -1 for opposite charge (K+ K-), +1 for same charge, 0 for any

# Resonance definitions
ChiCTo2Ka.states = cms.VPSet(
    cms.PSet(
        name = cms.string("ChiC0"),
        pdgId = cms.int32(10441),
        mass = cms.double(3.4147),
        massWindow = cms.double(0.080)
    ),
    cms.PSet(
        name = cms.string("ChiC2"),
        pdgId = cms.int32(445),
        mass = cms.double(3.5562),
        massWindow = cms.double(0.080)
    )
)
