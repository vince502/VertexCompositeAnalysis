"""
CMSSW_16_1_X MiniAOD configuration for conversion-based chi_c -> J/psi + gamma
producer and analyzer chains in VertexCompositeAnalysis.

Uses Era Run3_pp_on_PbPb_2026 (with automatic fallback to Run3_pp_on_PbPb_2025)
and GlobalTag auto:phase1_2026_realistic_hi by default, compatible with both
OniaTreeSubmodule (HiSkim/HiOnia2MuMu + HiAnalysis/HiOnia) and standalone CMSSW.
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

try:
    from Configuration.Eras.Era_Run3_pp_on_PbPb_2026_cff import (
        Run3_pp_on_PbPb_2026 as _DefaultHIEra,
    )
except ImportError:
    from Configuration.Eras.Era_Run3_pp_on_PbPb_2025_cff import (
        Run3_pp_on_PbPb_2025 as _DefaultHIEra,
    )

try:
    from Configuration.Eras.Era_Run3_2026_cff import Run3_2026 as _DefaultPPEra
except ImportError:
    from Configuration.Eras.Era_Run3_2025_cff import Run3_2025 as _DefaultPPEra

options = VarParsing("analysis")
options.register(
    "globalTag",
    "auto:phase1_2026_realistic_hi",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Conditions GlobalTag (161X / auto:phase1_2026_realistic_hi)",
)
options.register(
    "isMC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Enable MC truth branches and ChiCGenTruthNtuplizer",
)
options.register(
    "isPP",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use Run3 pp era instead of Run3_pp_on_PbPb HI era",
)
options.register(
    "doOneProng",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run single-leg conversion hypothesis producer and branches",
)
options.register(
    "doFlatNtuple",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run flat ChiCNtupleProducer alongside ChiCCovarianceNtuplizer",
)
options.register(
    "doHiOnia",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run companion HiOniaAnalyzer when OniaTreeSubmodule is installed",
)
options.outputFile = "chic_conversion_ntuple.root"
options.maxEvents = 100
options.parseArguments()

era = _DefaultPPEra if options.isPP else _DefaultHIEra
process = cms.Process("CHICCONVERSION", era)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, "")

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents),
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
)

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(options.outputFile),
)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCCombinedMiniAOD_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCEventInfoNtuplizer_cfi")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCNtupleProducer_cfi")

if not options.isPP:
    process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
    process.GlobalTag.toGet.append(process.chiCOfficialCentralityTable)

process.chiCCombinedCovarianceNtuplizer.isMC = cms.bool(options.isMC)
process.chiCCombinedCovarianceNtuplizer.doOneProng = cms.bool(options.doOneProng)
process.chiCEventInfoNtuplizer.isMC = cms.bool(options.isMC)
process.chiCGenTruthNtuplizer.doOneProng = cms.bool(options.doOneProng)

process.chiCNtupleProducer.candidateSrc = cms.InputTag("chiCConventionalPhotonStates")
process.chiCNtupleProducer.refitSrc = cms.InputTag("chiCConventionalFitter1S", "y1S")
process.chiCNtupleProducer.jpsiSrc = cms.InputTag("chiCSharedOniaJpsis")
process.chiCNtupleProducer.isCentrality = cms.untracked.bool(not options.isPP)
process.chiCNtupleProducer.centralityBinLabel = cms.InputTag("chiCOfficialCentralityBin", "HFtowers")
process.chiCNtupleProducer.isMC = cms.untracked.bool(options.isMC)
process.chiCNtupleProducer.storeRecoJpsi = cms.untracked.bool(True)

# Event-level info runs on all input events before any J/psi filter.
process.eventInfoPath = cms.Path(process.chiCEventInfoNtuplizer)

recoSeq = (
    process.chiCCombinedMiniAODSequence
    if options.doOneProng
    else process.chiCTwoProngMiniAODSequence
)

if not options.isPP:
    recoSeq = process.chiCOfficialCentralityBin * recoSeq

analyzerModules = [
    process.chiCCombinedCovarianceNtuplizer,
    process.oniaPhotonOfflineNtuplizer,
]

if options.doFlatNtuple:
    analyzerModules.append(process.chiCNtupleProducer)

if options.isMC:
    analyzerModules.append(process.chiCGenTruthNtuplizer)

if options.doHiOnia and getattr(process, "chiCSharedHiOnia", None) is not None:
    process.chiCSharedHiOnia.isHI = cms.untracked.bool(not options.isPP)
    process.chiCSharedHiOnia.isMC = cms.untracked.bool(options.isMC)
    analyzerModules.append(process.chiCSharedHiOnia)

analyzerSeq = analyzerModules[0]
for mod in analyzerModules[1:]:
    analyzerSeq = analyzerSeq * mod

process.chiCConversionPath = cms.Path(recoSeq * analyzerSeq)

process.schedule = cms.Schedule(
    process.eventInfoPath,
    process.chiCConversionPath,
)
