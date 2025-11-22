import FWCore.ParameterSet.Config as cms

DiMuonFromTracksCandidates = cms.EDProducer(
    "DiMuonFromTracksProducer",
    trackRecoAlgorithm = cms.InputTag("generalTracks"),
    vertexRecoAlgorithm = cms.InputTag("offlinePrimaryVertices"),
    tkChi2Cut = cms.double(10.0),
    tkNhitsCut = cms.int32(5),
    tkPtCut = cms.double(0.1),
    tkEtaCut = cms.double(2.4),
    tkDCACut = cms.double(1.0),
    mllCutMin = cms.double(2.5),  # J/ψ mass window
    mllCutMax = cms.double(3.5),
    jpsiMassCut = cms.double(0.15),  # ±150 MeV around J/ψ mass
    vtxChi2Cut = cms.double(10.0),
    dauTransImpactSigCut = cms.double(0.0),
    dauLongImpactSigCut = cms.double(0.0),
    doVertexFit = cms.bool(True),
    vertexFitter = cms.string("KalmanVertexFitter")
)
