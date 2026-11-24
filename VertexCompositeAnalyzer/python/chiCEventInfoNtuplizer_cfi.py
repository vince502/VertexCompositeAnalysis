import FWCore.ParameterSet.Config as cms

chiCEventInfoNtuplizer = cms.EDAnalyzer(
    "ChiCEventInfoNtuplizer",
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    beamSpot=cms.InputTag("offlineBeamSpot"),
    rho=cms.InputTag("fixedGridRhoFastjetAll"),
    hfRecHits=cms.InputTag("slimmedHcalRecHits", "reducedHcalRecHits"),
    centrality=cms.InputTag("hiCentrality"),
    pileupInfo=cms.InputTag("slimmedAddPileupInfo"),
    hltResults=cms.InputTag("TriggerResults", "", "HLT"),
    eventFilterResults=cms.InputTag("TriggerResults", "", "RECO"),
    trackMinPt=cms.double(0.4),
    trackMaxAbsEta=cms.double(2.4),
    trackMaxRelPtError=cms.double(0.1),
    trackMaxDzSignificance=cms.double(3.0),
    trackMaxDxySignificance=cms.double(3.0),
    hltPathPatterns=cms.vstring(
        "HLT_PPRefL*Mu*",
        "Dataset_PPRefDoubleMuon*",
        "Dataset_PPRefSingleMuon*",
    ),
    eventFilterPatterns=cms.vstring("Flag_*"),
    isMC=cms.bool(False),
)
