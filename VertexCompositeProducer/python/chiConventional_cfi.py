import FWCore.ParameterSet.Config as cms

chiOniaStates = cms.EDProducer(
    "ChiCOniaStateProducer",
    jpsis=cms.InputTag("onia2MuMuPatGlbGlb"),
    primaryVertices=cms.InputTag("unpackedTracksAndVertices"),
    useEmbeddedPVIndex=cms.bool(True),
    useOutputIndices=cms.bool(False),
    massRanges=cms.vdouble(2.6, 3.4),
    oniaPDGs=cms.vint32(443),
    nominalMasses=cms.vdouble(3.0969),
    oppositeSign=cms.bool(True),
)

chiConventionalCandidates = cms.EDProducer(
    "PoniaOniaPhotonProducer",
    dimuons=cms.InputTag("chiOniaStates"),
    conversions=cms.InputTag("oniaPhotonCandidates", "conversions"),
    pi0OnlineSwitch=cms.bool(False),
    deltaMass=cms.vdouble(0.0, 1.5),
    dzmax=cms.double(1.0),
    triggerMatch=cms.bool(False),
    minPhotonPt=cms.double(0.0),
    applyPairCuts=cms.bool(True),
    requireOniaVertex=cms.bool(True),
)

chiPhotonStates = cms.EDProducer(
    "ChiCConventionalPhotonStateProducer",
    candidates=cms.InputTag("chiConventionalCandidates"),
)

chiKinematicFit = cms.EDProducer(
    "OniaPhotonKinematicFit",
    chi_cand=cms.InputTag("chiPhotonStates"),
    upsilon_mass=cms.double(3.0969),
    product_name=cms.string("fit"),
    oniaPDGs=cms.vint32(443),
    nominalMasses=cms.vdouble(3.0969),
)
