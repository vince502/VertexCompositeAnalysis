import FWCore.ParameterSet.Config as cms

# Chi_c ntuple producer
# Reads chi_c candidates from OniaPhotonProducer / PoniaOniaPhotonProducer and writes to TTree
# Includes MC generator-level chi_c -> J/psi + gamma information
# Can also store raw reco::Conversion, PAT conversion, or Reco J/psi collections
chiCNtupleProducer = cms.EDAnalyzer(
    "ChiCNtupleProducer",
    candidateSrc=cms.InputTag("chiCandidatesMuMu"),
    pvSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
    # Centrality (for HI)
    isCentrality=cms.untracked.bool(False),
    centralityBinLabel=cms.InputTag("centralityBin", "HFtowers"),
    centralitySrc=cms.InputTag("hiCentrality"),
    genChiPDGs=cms.vint32(10441, 20443, 445, 10551, 20553, 555),
    genOniaPDGs=cms.vint32(443, 553),
    # MC generator info
    isMC=cms.untracked.bool(False),
    genParticles=cms.InputTag("prunedGenParticles"),
    # Reco conversions (photon -> e+e-)
    storeConversions=cms.untracked.bool(False),
    storePATConversions=cms.untracked.bool(False),
    patConversions=cms.InputTag("oniaPhotonCandidates", "conversions"),
    refitSrc=cms.InputTag(""),
    conversions=cms.InputTag("reducedEgamma", "reducedConversions"),
    # Reco J/psi collection
    storeRecoJpsi=cms.untracked.bool(False),
    jpsiSrc=cms.InputTag("onia2MuMuPAT"),
    # Tree name
    treeName=cms.untracked.string("chiTree"),
)

# Chi_c from mu+mu-
chiCMuMuNtupleProducer = chiCNtupleProducer.clone(
    candidateSrc=cms.InputTag("chiCandidatesMuMu"),
)

# Chi_c from e+e-
chiCEENtupleProducer = chiCNtupleProducer.clone(
    candidateSrc=cms.InputTag("chiCandidatesEE"),
)
