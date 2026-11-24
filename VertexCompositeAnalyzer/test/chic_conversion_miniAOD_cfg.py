"""
CMSSW_15_0_X MiniAOD configuration for conversion-based chi_c / chi_b -> Onia + gamma
reconstruction and ntuplization in 2025 OO, NeNe, pO (Run3_2025_OXY_cff) and
2025 PbPb (Run3_pp_on_PbPb_2025_cff) workflows.

Compatible with standalone VertexCompositeAnalysis and co-scheduled OniaTreeSubmodule
(HiSkim/HiOnia2MuMu + HiAnalysis/HiOnia on CMSSW_15_0_X).
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

try:
    from Configuration.Eras.Era_Run3_2025_OXY_cff import Run3_2025_OXY
except ImportError:
    from Configuration.StandardSequences.Eras import eras
    Run3_2025_OXY = eras.Run3_2025_OXY

try:
    from Configuration.Eras.Era_Run3_pp_on_PbPb_2025_cff import Run3_pp_on_PbPb_2025
except ImportError:
    try:
        from Configuration.Eras.Era_Run3_pp_on_PbPb_2024_cff import (
            Run3_pp_on_PbPb_2024 as Run3_pp_on_PbPb_2025,
        )
    except ImportError:
        from Configuration.Eras.Era_Run3_pp_on_PbPb_cff import (
            Run3_pp_on_PbPb as Run3_pp_on_PbPb_2025,
        )

options = VarParsing("analysis")
options.register(
    "collisionSystem",
    "OO",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Collision system: 'OO', 'NeNe', 'pO' (Run3_2025_OXY) or 'PbPb' (Run3_pp_on_PbPb_2025)",
)
options.register(
    "isMC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run MC truth matching and generator ntuplizers",
)
options.register(
    "doOneProng",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Enable single-leg conversion recovery branch alongside conventional 2-prong conversions",
)
options.register(
    "doFlatChiTree",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Schedule flat ChiCNtupleProducer alongside ChiCCovarianceNtuplizer",
)
options.register(
    "doHiOnia",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Schedule synchronized HiOniaAnalyzer (requires OniaTreeSubmodule HiAnalysis/HiOnia)",
)
options.register(
    "globalTag",
    "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Explicit GlobalTag override (defaults to 150X data/MC GT for selected collisionSystem)",
)
options.register(
    "ntupleOutputFile",
    "chic_conversion_150X_ntuple.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output ROOT file name for TFileService",
)
options.maxEvents = 100
options.parseArguments()

_system = options.collisionSystem.strip()
if _system in ("OO", "NeNe", "pO"):
    _era = Run3_2025_OXY
    _default_gt = "150X_mcRun3_2025_OO_realistic_v4" if options.isMC else "150X_dataRun3_Prompt_v1"
    _hlt_patterns = [
        "HLT_Oxy*Mu*",
        "HLT_MinimumBias*",
        "Dataset_Ion*",
    ]
elif _system == "PbPb":
    _era = Run3_pp_on_PbPb_2025
    _default_gt = "150X_mcRun3_2025_OO_realistic_v4" if options.isMC else "150X_dataRun3_Prompt_v1"
    _hlt_patterns = [
        "HLT_HIL*Mu*",
        "HLT_HIMinimumBias*",
        "Dataset_HIPhysicsRawPrime*",
    ]
else:
    raise ValueError(
        f"Unsupported collisionSystem='{_system}'. Choose from 'OO', 'NeNe', 'pO', or 'PbPb'."
    )

_resolved_gt = options.globalTag if options.globalTag else _default_gt

process = cms.Process("CHIC150X", _era)
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCCombinedMiniAOD_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCEventInfoNtuplizer_cfi")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCNtupleProducer_cfi")

process.GlobalTag = GlobalTag(process.GlobalTag, _resolved_gt, "")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
)
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
)
process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(options.ntupleOutputFile),
)

# For 2025 light-ion runs (OO/NeNe/pO), relax single-muon pT threshold to retain
# forward soft/tracker muons down to 1.5 GeV as in Run 3 OO dimuon workflows.
if _system in ("OO", "NeNe", "pO"):
    process.chiCMinPtSlimmedMuons.cut = cms.string("pt > 1.5")
    process.chiCSharedOniaSelectedMuons.cut = cms.string(
        "isTrackerMuon"
        ' && muonID("TMOneStationTight")'
        " && abs(innerTrack.dxy) < 0.3"
        " && abs(innerTrack.dz) < 20."
        " && innerTrack.hitPattern.trackerLayersWithMeasurement > 5"
        " && innerTrack.hitPattern.pixelLayersWithMeasurement > 0"
        " && abs(eta) < 2.4 && pt > 1.5"
    )
    process.chiCSharedOniaJpsis.dimuonSelection = cms.string(
        "2.6 < mass && mass < 3.5 && abs(y) < 2.4"
        " && charge == 0 && userFloat('vProb') > 0.01"
    )

process.chiCCombinedCovarianceNtuplizer.doOneProng = cms.bool(options.doOneProng)
process.chiCCombinedCovarianceNtuplizer.isMC = cms.bool(options.isMC)

process.chiCEventInfoNtuplizer.isMC = cms.bool(options.isMC)
process.chiCEventInfoNtuplizer.hltPathPatterns = cms.vstring(*_hlt_patterns)

process.offline = process.oniaPhotonOfflineNtuplizer.clone()

process.chiTree = process.chiCNtupleProducer.clone(
    candidateSrc=cms.InputTag("chiCConventionalPhotonStates"),
    pvSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
    refitSrc=cms.InputTag("chiCConventionalFitter1S", "y1S"),
    isCentrality=cms.untracked.bool(_system == "PbPb"),
    centralityBinLabel=cms.InputTag("chiCOfficialCentralityBin", "HFtowers"),
    centralitySrc=cms.InputTag("hiCentrality"),
    isMC=cms.untracked.bool(options.isMC),
    storePATConversions=cms.untracked.bool(False),
    storeRecoJpsi=cms.untracked.bool(True),
    jpsiSrc=cms.InputTag("chiCSharedOniaJpsis"),
)

_reco_steps = [
    process.chiCMinPtSlimmedMuonFilterSequence,
    process.chiCTwoProngMiniAODJpsiFilteredSequence,
]
if _system == "PbPb" or options.doHiOnia:
    _reco_steps.insert(0, process.chiCOfficialCentralityBin)
if options.doOneProng:
    _reco_steps.append(process.chiCOneProngMiniAODSequence)

_ana_steps = [
    process.chiCCombinedCovarianceNtuplizer,
    process.chiCEventInfoNtuplizer,
    process.offline,
]
if options.doFlatChiTree:
    _ana_steps.append(process.chiTree)
if options.isMC:
    process.chiCGenTruthNtuplizer.doOneProng = cms.bool(options.doOneProng)
    _ana_steps.append(process.chiCGenTruthNtuplizer)
if options.doHiOnia and hasattr(process, "chiCSharedHiOnia"):
    process.chiCSharedHiOnia.isMC = cms.untracked.bool(options.isMC)
    _ana_steps.append(process.chiCSharedHiOnia)

_full_seq = _reco_steps[0]
for _step in _reco_steps[1:] + _ana_steps:
    _full_seq = _full_seq * _step

process.p = cms.Path(_full_seq)
