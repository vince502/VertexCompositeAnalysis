import FWCore.ParameterSet.Config as cms

trackCountFilter = cms.EDFilter("TrackCountFilter",
    trackCollection = cms.InputTag("generalTracks"),
    minTrackPt = cms.double(0.4),
    maxTrackEta = cms.double(2.4),
    maxTrackNormalizedChi2 = cms.double(10.0),
    minTrackNHits = cms.int32(6),
    minTrackNPix = cms.int32(0),  # 0 = no cut (can be overridden)
    minNSelectedTracks = cms.int32(4),
    maxNSelectedTracks = cms.int32(999999)
)
