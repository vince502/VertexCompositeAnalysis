"""CMSSW_13_2_X MiniAOD configuration for conversion-based chi_c -> J/psi + gamma production & ntuplization."""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
from Configuration.Eras.Modifier_run3_upc_cff import run3_upc

options = VarParsing("analysis")
options.register(
    "isMC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on MiniAODSIM (MC) instead of Data MiniAOD",
)
options.register(
    "isHI",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use 2023 PbPb heavy-ion era and centrality settings",
)
options.register(
    "isUPC",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply Run 3 UPC modifier (Run3_2023 + run3_upc)",
)
options.register(
    "doOneProng",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Enable single-leg conversion hypothesis branch alongside two-prong conversions",
)
options.register(
    "doFlatNtuple",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run flat ChiCNtupleProducer alongside ChiCCovarianceNtuplizer",
)
options.register(
    "useStandaloneJpsi",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Build dimuons with standalone ChiCJpsiProducer instead of Onia2MuMuPAT + DiMuonFilter",
)
options.register(
    "globalTag",
    "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Override GlobalTag (defaults to 132X_dataRun3_Prompt_v7 or auto:phase1_2023_realistic_hi)",
)
options.maxEvents = 100
options.outputFile = "chic_conversion_ntuple_132X.root"
if not options.inputFiles:
    options.inputFiles = [
        "file:/eos/cms/store/group/phys_heavyions/dileptons/Data2023/MINIAOD/HIPhysicsRawPrime0/Run375064/7ed5766f-6b1d-415e-8916-e62825a6347f.root"
    ]
options.parseArguments()

if options.isUPC:
    process = cms.Process("CHICCONV", Run3_2023, run3_upc)
elif options.isHI:
    process = cms.Process("CHICCONV", Run3_pp_on_PbPb_2023)
else:
    process = cms.Process("CHICCONV", Run3_2023)

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

resolvedGT = options.globalTag
if not resolvedGT:
    resolvedGT = "auto:phase1_2023_realistic_hi" if options.isMC else "132X_dataRun3_Prompt_v7"
process.GlobalTag = GlobalTag(process.GlobalTag, resolvedGT, "")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
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

if options.isHI:
    process.GlobalTag.toGet.append(process.chiCOfficialCentralityTable)
    process.chiCEventInfoNtuplizer.hltPathPatterns = cms.vstring(
        "HLT_HIL*Mu*",
        "HLT_HIUPC*Mu*",
        "HLT_HIMinimumBias*",
    )

process.chiCEventInfoNtuplizer.isMC = cms.bool(options.isMC)
process.chiCCombinedCovarianceNtuplizer.isMC = cms.bool(options.isMC)
process.chiCCombinedCovarianceNtuplizer.doOneProng = cms.bool(options.doOneProng)
process.chiCGenTruthNtuplizer.doOneProng = cms.bool(options.doOneProng)

if options.useStandaloneJpsi:
    process.chiCSharedOniaJpsis = process.chiCJpsis.clone(
        muons=cms.InputTag("chiCSharedOniaSelectedMuons"),
    )
    process.chiCSharedOniaMiniAODJpsiFilteredSequence = cms.Sequence(
        process.chiCSharedOniaSelectedMuons
        * process.chiCSharedOniaJpsis
        * process.chiCSharedOniaJpsiCounter
        * process.chiCSharedOniaStates
    )

recoSeq = cms.Sequence()
if options.isHI:
    recoSeq += process.chiCOfficialCentralityBin
recoSeq += process.chiCMinPtSlimmedMuonFilterSequence
recoSeq += process.chiCSharedOniaMiniAODJpsiFilteredSequence
recoSeq += process.chiCConventionalMiniAODSequence
if options.doOneProng:
    recoSeq += process.chiCOneProngMiniAODSequence

anaSeq = cms.Sequence(
    process.chiCCombinedCovarianceNtuplizer
    * process.chiCEventInfoNtuplizer
    * process.oniaPhotonOfflineNtuplizer
)

if options.isMC:
    anaSeq += process.chiCGenTruthNtuplizer

if options.doFlatNtuple:
    process.chiCNtupleProducer.candidateSrc = cms.InputTag("chiCConventionalPhotonStates")
    process.chiCNtupleProducer.refitSrc = cms.InputTag("chiCConventionalFitter1S", "y1S")
    process.chiCNtupleProducer.isCentrality = cms.untracked.bool(options.isHI)
    process.chiCNtupleProducer.centralityBinLabel = cms.InputTag("chiCOfficialCentralityBin", "HFtowers")
    process.chiCNtupleProducer.isMC = cms.untracked.bool(options.isMC)
    anaSeq += process.chiCNtupleProducer

if hasattr(process, "chiCSharedHiOnia") and options.isHI and not options.useStandaloneJpsi:
    process.chiCSharedHiOnia.isMC = cms.untracked.bool(options.isMC)
    anaSeq += process.chiCSharedHiOnia

process.p = cms.Path(recoSeq * anaSeq)
