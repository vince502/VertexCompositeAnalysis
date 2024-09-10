import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_36_patch1/src/Closure/HIN-pPb816Spring16GS-00122.root',
#'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/CMSSW_8_0_36_patch1/src/Closure/HIN-pPb816Spring16GS-00122_Soft.root',
'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/EvtGenMod/CMSSW_8_0_30/src/Closures/HIN-pPb816Spring16GS-00122.root',
)
)

# =============== Other Statements =====================
# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '80X_dataRun2_v19'

# =============== Import Sequences =====================

process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.simplegen_cfi')
process.analysis = cms.Path(process.genAna);
process.schedule = cms.Schedule(process.analysis)

#process.output = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM", "keep *_*_generalTracks_*"),
#    fileName = cms.untracked.string('output.root'),
#)
#
#process.outputPath = cms.EndPath(process.output)
#process.schedule.append(process.outputPath)
