import FWCore.ParameterSet.Config as cms

ChiCTo4Pi = cms.EDProducer("ChiCFourTrackProducer")

# Track selection cuts
ChiCTo4Pi.trackCollection = cms.InputTag("generalTracks")
ChiCTo4Pi.daughterMasses = cms.vdouble(0.13957018, 0.13957018, 0.13957018, 0.13957018)
ChiCTo4Pi.minTrackPt = cms.double(0.4)
ChiCTo4Pi.maxTrackEta = cms.double(2.4)
ChiCTo4Pi.maxTrackNormalizedChi2 = cms.double(10.0)
ChiCTo4Pi.minTrackNHits = cms.int32(6)

# Candidate-level cuts
ChiCTo4Pi.minCandidatePt = cms.double(3.5)
ChiCTo4Pi.minAcoplanarity = cms.double(0.6)
ChiCTo4Pi.maxSphericity = cms.double(0.35)
ChiCTo4Pi.applyMassWindow = cms.bool(True)

# Resonance definitions
ChiCTo4Pi.states = cms.VPSet(
    cms.PSet(
        name = cms.string("ChiC0"),
        pdgId = cms.int32(10441),
        mass = cms.double(3.4147),
        massWindow = cms.double(0.025)
    ),
    cms.PSet(
        name = cms.string("ChiC2"),
        pdgId = cms.int32(445),
        mass = cms.double(3.5562),
        massWindow = cms.double(0.025)
    )
)
