import FWCore.ParameterSet.Config as cms

chiCCovarianceNtuplizer = cms.EDAnalyzer(
    "ChiCCovarianceNtuplizer",
    oniaStates=cms.InputTag("chiCSharedOniaStates"),
    oniaMuons=cms.InputTag("chiCSharedOniaSelectedMuons"),
    conventionalCandidates=cms.InputTag("chiCConventionalCandidates"),
    oneProngCandidates=cms.InputTag("chiCSharedOniaOneProngCandidates"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    trackMinPt=cms.double(0.4),
    trackMaxAbsEta=cms.double(2.4),
    trackMaxRelPtError=cms.double(0.1),
    trackMaxDzSignificance=cms.double(3.0),
    trackMaxDxySignificance=cms.double(3.0),
    doOneProng=cms.bool(True),
    isMC=cms.bool(False),
)

chiCGenTruthNtuplizer = cms.EDAnalyzer(
    "ChiCGenTruthNtuplizer",
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    oniaStates=cms.InputTag("chiCSharedOniaStates"),
    conventionalCandidates=cms.InputTag("chiCConventionalPhotonStates"),
    oneProngCandidates=cms.InputTag("chiCSharedOniaOneProngCandidates"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    packedToGenAssociation=cms.InputTag("packedPFCandidateToGenAssociation"),
    lostToGenAssociation=cms.InputTag("lostTracksToGenAssociation"),
    doOneProng=cms.bool(True),
)

oniaPhotonOfflineNtuplizer = cms.EDAnalyzer(
    "OniaPhotonOfflineNtuplizer",
    onias=cms.InputTag("chiCSharedOniaStates"),
    chis=cms.InputTag("chiCConventionalPhotonStates"),
    muons=cms.InputTag("slimmedMuons"),
    vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    trackMinPt=cms.double(0.4),
    trackMaxAbsEta=cms.double(2.4),
    trackMaxRelPtError=cms.double(0.1),
    trackMaxDzSignificance=cms.double(3.0),
    trackMaxDxySignificance=cms.double(3.0),
    triggerObjects=cms.InputTag("slimmedPatTrigger"),
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
)

chiCOneProngInputAudit = cms.EDAnalyzer(
    "ChiCOneProngInputAudit",
    muons=cms.InputTag("slimmedMuons"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot=cms.InputTag("offlineBeamSpot"),
)
