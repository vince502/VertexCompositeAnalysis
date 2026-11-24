"""
CMSSW_14_1_X MiniAOD configuration template for chi_c -> J/psi + gamma
conversion-based photon producer and analyzer workflows (2024 ppRef & 2024 PbPb).

Examples:
  # 2024 ppRef PromptReco Data (default):
  cmsRun chic_conversion_miniAOD_cfg.py inputFiles=file:input_miniAOD.root

  # 2024 ppRef MC (with one-prong & generator truth):
  cmsRun chic_conversion_miniAOD_cfg.py isMC=True doOneProng=True inputFiles=file:mc_miniAOD.root

  # 2024 PbPb PromptReco Data:
  cmsRun chic_conversion_miniAOD_cfg.py collisionSystem=PbPb2024 inputFiles=file:pbpb_miniAOD.root
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

options = VarParsing("analysis")
options.register(
    "collisionSystem",
    "ppRef2024",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Collision system: 'ppRef2024' or 'PbPb2024'",
)
options.register(
    "isMC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run MC truth matching and pileup summary branches",
)
options.register(
    "doOneProng",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run single-leg conversion recovery branch alongside two-prong conversions",
)
options.register(
    "doFlatNtuple",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Also run flat ChiCNtupleProducer tree alongside covariance/offline trees",
)
options.register(
    "doHiOnia",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Schedule companion HiOniaAnalyzer tree when OniaTreeSubmodule is installed",
)
options.register(
    "broadOnia",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use broad OniaMuonSelector (pt > 2 GeV, offline IDs) instead of standard soft-muon preselection",
)
options.register(
    "globalTag",
    "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Override GlobalTag (defaults to 141X data/MC tags per collisionSystem)",
)
options.maxEvents = 100
options.outputFile = "chic_conversion_ntuple.root"
options.parseArguments()

is_pbpb = options.collisionSystem.lower() in ("pbpb2024", "pbpb", "hi")

if is_pbpb:
    try:
        from Configuration.Eras.Era_Run3_pp_on_PbPb_2024_cff import (
            Run3_pp_on_PbPb_2024 as selectedEra,
        )
    except ImportError:
        from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import (
            Run3_pp_on_PbPb_2023 as selectedEra,
        )
    default_gt = (
        "auto:phase1_2024_realistic_hi" if options.isMC else "141X_dataRun3_Prompt_v3"
    )
else:
    from Configuration.Eras.Era_Run3_2024_ppRef_cff import Run3_2024_ppRef as selectedEra

    default_gt = (
        "141X_mcRun3_2024_realistic_ppRef5TeV_v7"
        if options.isMC
        else "141X_dataRun3_Prompt_v3"
    )

resolved_gt = options.globalTag if options.globalTag else default_gt

process = cms.Process("CHICCONV", selectedEra)
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, resolved_gt, "")

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents)
)
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
)
process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(options.outputFile),
)

if is_pbpb:
    process.load(
        "VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCCombinedMiniAOD_cff"
    )
    process.load(
        "VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCEventInfoNtuplizer_cfi"
    )
    process.GlobalTag.toGet.extend([process.chiCOfficialCentralityTable])
    process.chiCEventInfoNtuplizer.isMC = cms.bool(options.isMC)
    process.chiCEventInfoNtuplizer.hltPathPatterns = cms.vstring(
        "HLT_HIL*Mu*",
        "Dataset_HIPhysicsRawPrime*",
    )
    process.chiCCombinedCovarianceNtuplizer.isMC = cms.bool(options.isMC)
    process.chiCCombinedCovarianceNtuplizer.doOneProng = cms.bool(options.doOneProng)
    reco_seq = cms.Sequence(
        process.chiCOfficialCentralityBin
        * process.chiCMinPtSlimmedMuonFilterSequence
        * process.chiCTwoProngMiniAODJpsiFilteredSequence
    )
    if options.doOneProng:
        reco_seq *= process.chiCOneProngMiniAODSequence
    cov_ntuplizer = process.chiCCombinedCovarianceNtuplizer
    evt_ntuplizer = process.chiCEventInfoNtuplizer
else:
    process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCPPRefNtuple_cff")
    if options.broadOnia:
        old_mu_sel = process.chiCSharedOniaSelectedMuons
        process.chiCSharedOniaSelectedMuons = cms.EDProducer(
            "OniaMuonSelector",
            src=cms.InputTag("slimmedMuons"),
            minPt=cms.double(2.0),
        )
        process.chiCSharedOniaMiniAODSequence.replace(
            old_mu_sel, process.chiCSharedOniaSelectedMuons
        )
        process.chiCSharedOniaMiniAODJpsiFilteredSequence.replace(
            old_mu_sel, process.chiCSharedOniaSelectedMuons
        )
        process.chiCSharedOniaJpsis.dimuonSelection = cms.string(
            "2.5 < mass && mass < 3.6 && abs(y) < 2.5 && userFloat('vProb') > 0."
        )
        process.chiCSharedOniaStates.useEmbeddedPVIndex = cms.bool(False)
        process.chiCSharedOniaStates.useOutputIndices = cms.bool(True)
        process.chiCConventionalCandidates.applyPairCuts = cms.bool(False)
        process.chiCConventionalCandidates.requireOniaVertex = cms.bool(True)
    if options.isMC:
        process.chiCPPRefMCCovarianceNtuplizer.doOneProng = cms.bool(options.doOneProng)
        reco_seq = cms.Sequence(
            process.chiCMinPtSlimmedMuonFilterSequence
            * process.chiCTwoProngMiniAODJpsiFilteredSequence
        )
        if options.doOneProng:
            reco_seq *= process.chiCOneProngMiniAODSequence
        cov_ntuplizer = process.chiCPPRefMCCovarianceNtuplizer
        evt_ntuplizer = process.chiCPPRefMCEventInfoNtuplizer
    else:
        process.chiCPPRefDataCovarianceNtuplizer.doOneProng = cms.bool(
            options.doOneProng
        )
        reco_seq = cms.Sequence(process.chiCPPRefDataSequence)
        if options.doOneProng:
            reco_seq *= process.chiCOneProngMiniAODSequence
        cov_ntuplizer = process.chiCPPRefDataCovarianceNtuplizer
        evt_ntuplizer = process.chiCPPRefDataEventInfoNtuplizer

from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCCovarianceNtuplizer_cfi import (
    oniaPhotonOfflineNtuplizer,
)

process.offline = oniaPhotonOfflineNtuplizer.clone()

analysis_seq = cms.Sequence(reco_seq * cov_ntuplizer * evt_ntuplizer * process.offline)

if options.isMC:
    process.chiCGenTruthNtuplizer.doOneProng = cms.bool(options.doOneProng)
    analysis_seq *= process.chiCGenTruthNtuplizer

if options.doFlatNtuple:
    from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCNtupleProducer_cfi import (
        chiCNtupleProducer,
    )

    process.chiFlatNtuple = chiCNtupleProducer.clone(
        candidateSrc=cms.InputTag("chiCConventionalPhotonStates"),
        pvSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
        isCentrality=cms.untracked.bool(is_pbpb),
        centralityBinLabel=cms.InputTag("chiCOfficialCentralityBin", "HFtowers"),
        centralitySrc=cms.InputTag("hiCentrality"),
        isMC=cms.untracked.bool(options.isMC),
        refitSrc=cms.InputTag("chiCConventionalFitter1S", "y1S"),
        storeRecoJpsi=cms.untracked.bool(True),
        jpsiSrc=cms.InputTag("chiCSharedOniaStates"),
    )
    analysis_seq *= process.chiFlatNtuple

if options.doHiOnia and hasattr(process, "chiCSharedHiOnia"):
    process.chiCSharedHiOnia.isHI = cms.untracked.bool(is_pbpb)
    process.chiCSharedHiOnia.isMC = cms.untracked.bool(options.isMC)
    process.chiCSharedHiOnia.useEvtPlane = cms.untracked.bool(is_pbpb)
    if not is_pbpb:
        process.chiCSharedHiOnia.CentralitySrc = cms.InputTag("")
        process.chiCSharedHiOnia.CentralityBinSrc = cms.InputTag("")
    analysis_seq *= process.chiCSharedHiOnia

process.p = cms.Path(analysis_seq)
