import FWCore.ParameterSet.Config as cms

ChiC2To4K = cms.EDProducer("ChiCFourTrackProducer")

# Track selection cuts
ChiC2To4K.trackCollection = cms.InputTag("generalTracks")
ChiC2To4K.daughterMasses = cms.vdouble(0.493677, 0.493677, 0.493677, 0.493677)
ChiC2To4K.minTrackPt = cms.double(0.6)
ChiC2To4K.maxTrackEta = cms.double(2.4)
ChiC2To4K.maxTrackNormalizedChi2 = cms.double(10.0)
ChiC2To4K.minTrackNHits = cms.int32(6)

# Candidate-level cuts
ChiC2To4K.minCandidatePt = cms.double(1.5)
ChiC2To4K.minAcoplanarity = cms.double(0.0)
ChiC2To4K.maxSphericity = cms.double(1.5)
ChiC2To4K.maxCandidateAbsEta = cms.double(2.4)
ChiC2To4K.storeEventShape = cms.bool(False)
ChiC2To4K.applyMassWindow = cms.bool(True)

# Resonance definitions
ChiC2To4K.states = cms.VPSet(
    cms.PSet(
        name = cms.string("ChiC0"),
        pdgId = cms.int32(10441),
        mass = cms.double(3.4147),
        massWindow = cms.double(0.090)
    ),
    cms.PSet(
        name = cms.string("ChiC2"),
        pdgId = cms.int32(445),
        mass = cms.double(3.5562),
        massWindow = cms.double(0.090)
    )
)
