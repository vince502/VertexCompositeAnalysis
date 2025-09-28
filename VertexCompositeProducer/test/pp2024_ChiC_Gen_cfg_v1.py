import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('GENANA', eras.Run3_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(1),
)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        # Replace with GEN-SIM-RECO/MINIAODSIM files that contain prunedGenParticles
        'file:/afs/cern.ch/work/s/soohwan/private/Analysis/MC/ppRef2024/CMSSW_14_1_9/src/ChicTo4Pi.root'
    ),
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('141X_dataRun3_Express_v3')

from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCGenNtuplizer_cfi import ChiCGenNtuplizer as _ChiCGenNtuplizer

process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string('chic_gen_ntuple.root'),
)

process.GenChiCNtuplizer = _ChiCGenNtuplizer.clone(
    genParticles=cms.InputTag('genParticles'),
    primaryVertices=cms.InputTag('offlineSlimmedPrimaryVertices'),
    treeName=cms.untracked.string('ChiCGenNtuple'),
)

process.chicGenEndPath = cms.EndPath(process.GenChiCNtuplizer)

process.schedule = cms.Schedule(process.chicGenEndPath)
