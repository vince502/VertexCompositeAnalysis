import FWCore.ParameterSet.Config as cms

try:
    from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cfi import onia2MuMuPAT
except ImportError:
    from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi import onia2MuMuPAT

try:
    from HiAnalysis.HiOnia.hioniaanalyzer_cfi import hionia as _hionia_base
except ImportError:
    _hionia_base = None

try:
    from RecoHI.HiCentralityAlgos.CentralityBin_cfi import centralityBin as _centrality_bin_base
except ImportError:
    _centrality_bin_base = None

from VertexCompositeAnalysis.VertexCompositeProducer.chiCOneProng_cfi import (
    chiCConversionLegs,
    chiCConventionalPhotonStates as _chiCConventionalPhotonStatesBase,
    chiCJpsis,
    chiCOniaStates,
    chiCOneProngCandidates,
    chiCSingleLegGeometries,
)
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCCovarianceNtuplizer_cfi import (
    chiCCovarianceNtuplizer,
    chiCGenTruthNtuplizer,
    chiCOneProngInputAudit,
    oniaPhotonOfflineNtuplizer,
)

chiCOfficialCentralityTable = cms.PSet(
    record=cms.string("HeavyIonRcd"),
    tag=cms.string(
        "CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_Nominal"
    ),
    label=cms.untracked.string("HFtowers"),
)

if _centrality_bin_base is not None:
    chiCOfficialCentralityBin = _centrality_bin_base.clone(
        Centrality=cms.InputTag("hiCentrality"),
        centralityVariable=cms.string("HFtowers"),
        nonDefaultGlauberModel=cms.string(""),
    )
else:
    chiCOfficialCentralityBin = cms.EDProducer(
        "CentralityBinProducer",
        Centrality=cms.InputTag("hiCentrality"),
        centralityVariable=cms.string("HFtowers"),
        nonDefaultGlauberModel=cms.string(""),
        pPbRunFlip=cms.uint32(99999999),
        effAndAsymRun=cms.uint32(0),
    )

chiCMinPtSlimmedMuons = cms.EDFilter(
    "PATMuonSelector",
    src=cms.InputTag("slimmedMuons"),
    cut=cms.string("pt > 2.5"),
    filter=cms.bool(False),
)

chiCMinPtSlimmedMuonCounter = cms.EDFilter(
    "CandViewCountFilter",
    src=cms.InputTag("chiCMinPtSlimmedMuons"),
    minNumber=cms.uint32(2),
)

# This is the one shared J/psi collection for both photon branches.
chiCSharedOniaSelectedMuons = cms.EDFilter(
    "PATMuonSelector",
    src=cms.InputTag("slimmedMuons"),
    cut=cms.string(
        "isGlobalMuon && isTrackerMuon"
        ' && muonID("TMOneStationTight")'
        " && abs(innerTrack.dxy) < 0.3"
        " && abs(innerTrack.dz) < 20."
        " && innerTrack.hitPattern.trackerLayersWithMeasurement > 5"
        " && innerTrack.hitPattern.pixelLayersWithMeasurement > 0"
        ' && innerTrack.quality("highPurity")'
        " && abs(eta) < 2.4 && pt > 3.5"
    ),
    filter=cms.bool(False),
)

chiCSharedOnia2MuMu = onia2MuMuPAT.clone(
    muons=cms.InputTag("chiCSharedOniaSelectedMuons"),
    primaryVertexTag=cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpotTag=cms.InputTag("offlineBeamSpot"),
    higherPuritySelection=cms.string(""),
    lowerPuritySelection=cms.string(""),
    dimuonSelection=cms.string("2.7 < mass && mass < 3.5"),
)
if hasattr(chiCSharedOnia2MuMu, "addMCTruth"):
    chiCSharedOnia2MuMu.addMCTruth = cms.bool(False)

chiCSharedOniaJpsis = cms.EDProducer(
    "DiMuonFilter",
    OniaTag=cms.InputTag("chiCSharedOnia2MuMu"),
    singlemuonSelection=cms.string(""),
    dimuonSelection=cms.string(
        "2.7 < mass && mass < 3.5 && pt > 4. && abs(y) < 2.4"
        " && charge == 0 && userFloat('vProb') > 0.005"
    ),
    do_trigger_match=cms.bool(False),
    HLTFilters=cms.vstring(),
)

chiCSharedOniaJpsiCounter = cms.EDFilter(
    "CandViewCountFilter",
    src=cms.InputTag("chiCSharedOniaJpsis"),
    minNumber=cms.uint32(1),
)

# Explicit vertex-fitted J/psi states shared by both photon branches.
chiCSharedOniaStates = chiCOniaStates.clone(
    jpsis=cms.InputTag("chiCSharedOniaJpsis"),
)

# Conventional fully reconstructed gamma -> e+e- branch.
chiCConventionalCandidates = cms.EDProducer(
    "PoniaOniaPhotonProducer",
    conversions=cms.InputTag("oniaPhotonCandidates", "conversions"),
    dimuons=cms.InputTag("chiCSharedOniaStates"),
    pi0OnlineSwitch=cms.bool(False),
    deltaMass=cms.vdouble(0.0, 1.5),
    dzmax=cms.double(1.0),
    triggerMatch=cms.bool(False),
    minPhotonPt=cms.double(0.0),
    applyPairCuts=cms.bool(True),
    requireOniaVertex=cms.bool(True),
)

chiCConventionalPhotonStates = _chiCConventionalPhotonStatesBase.clone(
    candidates=cms.InputTag("chiCConventionalCandidates"),
)

chiCConventionalFitter1S = cms.EDProducer(
    "OniaPhotonKinematicFit",
    chi_cand=cms.InputTag("chiCConventionalPhotonStates"),
    upsilon_mass=cms.double(3.0969),
    product_name=cms.string("y1S"),
    oniaPDGs=cms.vint32(),
    nominalMasses=cms.vdouble(),
)

# One-prong branch. Its photon hypotheses consume the shared fitted J/psi
# states used by the conventional branch.
chiCSharedOniaConversionLegs = chiCConversionLegs.clone()
chiCSharedOniaSingleLegGeometries = chiCSingleLegGeometries.clone(
    jpsis=cms.InputTag("chiCSharedOniaStates"),
    legs=cms.InputTag("chiCSharedOniaConversionLegs"),
)
chiCSharedOniaOneProngCandidates = chiCOneProngCandidates.clone(
    geometries=cms.InputTag("chiCSharedOniaSingleLegGeometries"),
    jpsis=cms.InputTag("chiCSharedOniaStates"),
    legs=cms.InputTag("chiCSharedOniaConversionLegs"),
)

chiCCombinedCovarianceNtuplizer = chiCCovarianceNtuplizer.clone(
    oniaStates=cms.InputTag("chiCSharedOniaStates"),
    oniaMuons=cms.InputTag("chiCSharedOniaSelectedMuons"),
    conventionalCandidates=cms.InputTag("chiCConventionalPhotonStates"),
    oneProngCandidates=cms.InputTag("chiCSharedOniaOneProngCandidates"),
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    doOneProng=cms.bool(True),
    isMC=cms.bool(False),
)

# When OniaTreeSubmodule (HiAnalysis/HiOnia) is present, provide chiCSharedHiOnia.
if _hionia_base is not None:
    chiCSharedHiOnia = _hionia_base.clone(
        srcMuon=cms.InputTag("chiCSharedOniaSelectedMuons"),
        srcMuonNoTrig=cms.InputTag("chiCSharedOniaSelectedMuons"),
        srcDimuon=cms.InputTag("chiCSharedOniaJpsis"),
        primaryVertexTag=cms.InputTag("offlineSlimmedPrimaryVertices"),
        triggerResultsLabel=cms.InputTag("TriggerResults", "", "HLT"),
        stageL1Trigger=cms.uint32(2),
        CentralitySrc=cms.InputTag("hiCentrality"),
        CentralityBinSrc=cms.InputTag("chiCOfficialCentralityBin", "HFtowers"),
        EvtPlane=cms.InputTag("hiEvtPlaneFlat"),
        maxAbsZ=cms.double(1.0e6),
        muonSel=cms.string("All"),
        isHI=cms.untracked.bool(True),
        isPA=cms.untracked.bool(False),
        isMC=cms.untracked.bool(False),
        applyCuts=cms.bool(False),
        miniAODcut=cms.bool(False),
        checkTrigNames=cms.bool(False),
        SumETvariables=cms.bool(True),
        useEvtPlane=cms.untracked.bool(True),
        useGeTracks=cms.untracked.bool(False),
        fillRecoTracks=cms.bool(False),
        fillSingleMuons=cms.bool(True),
        AtLeastOneCand=cms.bool(False),
        mom4format=cms.string("vector"),
        dblTriggerPathNames=cms.vstring(),
        sglTriggerPathNames=cms.vstring(),
    )

chiCMinPtSlimmedMuonFilterSequence = cms.Sequence(
    chiCMinPtSlimmedMuons
    * chiCMinPtSlimmedMuonCounter
)

chiCSharedOniaMiniAODSequence = cms.Sequence(
    chiCSharedOniaSelectedMuons
    * chiCSharedOnia2MuMu
    * chiCSharedOniaJpsis
    * chiCSharedOniaStates
)

chiCSharedOniaMiniAODJpsiFilteredSequence = cms.Sequence(
    chiCSharedOniaSelectedMuons
    * chiCSharedOnia2MuMu
    * chiCSharedOniaJpsis
    * chiCSharedOniaJpsiCounter
    * chiCSharedOniaStates
)

chiCConventionalMiniAODSequence = cms.Sequence(
    chiCConventionalCandidates
    * chiCConventionalPhotonStates
    * chiCConventionalFitter1S
)

chiCOneProngMiniAODSequence = cms.Sequence(
    chiCSharedOniaConversionLegs
    * chiCSharedOniaSingleLegGeometries
    * chiCSharedOniaOneProngCandidates
)

chiCTwoProngMiniAODSequence = cms.Sequence(
    chiCSharedOniaMiniAODSequence
    * chiCConventionalMiniAODSequence
)

chiCTwoProngMiniAODJpsiFilteredSequence = cms.Sequence(
    chiCSharedOniaMiniAODJpsiFilteredSequence
    * chiCConventionalMiniAODSequence
)

chiCCombinedMiniAODSequence = cms.Sequence(
    chiCTwoProngMiniAODSequence
    * chiCOneProngMiniAODSequence
)
