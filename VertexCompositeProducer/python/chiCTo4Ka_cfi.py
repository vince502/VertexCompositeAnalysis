import FWCore.ParameterSet.Config as cms

ChiCTo4Ka = cms.EDProducer(
    "ChiCFourTrackProducer",
    trackCollection = cms.InputTag("generalTracks"),
    daughterMasses = cms.vdouble(0.493677, 0.493677, 0.493677, 0.493677),
    minTrackPt = cms.double(0.6),
    maxTrackEta = cms.double(2.4),
    maxTrackNormalizedChi2 = cms.double(10.0),
    minTrackNHits = cms.int32(6),
    minCandidatePt = cms.double(1.5),
    minAcoplanarity = cms.double(0.0),
    maxSphericity = cms.double(1.5),
    maxCandidateAbsEta = cms.double(2.4),
    storeEventShape = cms.bool(False),
    applyMassWindow = cms.bool(True),
    states = cms.VPSet(
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
)
