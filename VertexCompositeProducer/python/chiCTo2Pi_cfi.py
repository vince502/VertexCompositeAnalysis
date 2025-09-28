import FWCore.ParameterSet.Config as cms

ChiCTo2Pi = cms.EDProducer(
    "ChiCTrackPairProducer",
    trackCollection = cms.InputTag("generalTracks"),
    daughterMass = cms.double(0.13957018),
    minTrackPt = cms.double(0.5),
    maxTrackEta = cms.double(2.4),
    maxTrackNormalizedChi2 = cms.double(10.0),
    minTrackNHits = cms.int32(6),
    minPairPt = cms.double(1.0),
    applyMassWindow = cms.bool(True),
    states = cms.VPSet(
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
)
