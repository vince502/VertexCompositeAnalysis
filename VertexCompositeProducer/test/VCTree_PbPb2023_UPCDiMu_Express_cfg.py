import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run3_pp_on_PbPb)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
# process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True), SkipEvent = cms.untracked.vstring('ProductNotFound'))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/wangj/RECO2023/aod_PhysicsHIForward0_374345/reco_run374345_ls0100_streamPhysicsHIForward0_StorageManager.root'),
   inputCommands=cms.untracked.vstring('keep *', 'drop *_hiEvtPlane_*_*'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Express_v4')


## ##############################################################################################################################
## Variables Production #########################################################################################################

#* cent_seq: Add PbPb centrality
#* CentralityBin_cfi contains Ntrack information
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        # tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline"),
        # tag = cms.string("CentralityTable_HFtowers200_v0_express"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.centralityBin)

#* Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, False)

#* Muon selection
process.twoMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"), minNumber = cms.uint32(2))
process.goodMuon = cms.EDFilter("MuonSelector",
            src = cms.InputTag("muons"),
            cut = process.generalMuMuMassMin0Candidates.muonSelection,
            )
process.twoGoodMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodMuon"), minNumber = cms.uint32(2))
process.goodDimuon = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = process.generalMuMuMassMin0Candidates.candidateSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodMuon@+ goodMuon@-')
            )
process.oneGoodDimuon = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDimuon"), minNumber = cms.uint32(1))
process.dimuonEvtSel = cms.Sequence(process.twoMuons * process.goodMuon * process.twoGoodMuons * process.goodDimuon * process.oneGoodDimuon)

#* Add PbPb collision event selection
#* hfCoincFilter2Th4: 2 HF towers above 4 GeV
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.clusterCompatibilityFilter_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load("VertexCompositeAnalysis.VertexCompositeProducer.OfflinePrimaryVerticesRecovery_cfi")
# process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter)
process.colEvtSel = cms.Sequence(process.primaryVertexFilter)

#* Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
	process.colEvtSel *
    process.dimuonEvtSel *
	process.offlinePrimaryVerticesRecovery
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#* Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq)
process.dimurereco_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin0Candidates)


## Adding the VertexComposite tree ################################################################################################

# Add the VertexComposite ntuple
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")
process.dimucontana.selectEvents = cms.untracked.string("eventFilter_HM_step")
process.dimucontana.doMuonNtuple = cms.untracked.bool(True)
process.dimucontana.VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin0Candidates:DiMu")
# process.dimucontana_wrongsign.selectEvents = cms.untracked.string("eventFilter_HM_step")

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana.root'))
process.p = cms.EndPath(process.dimucontana)

#! Define the process schedule !!!!!!!!!!!!!!!!!!
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    # process.pcentandep_step,
    process.dimurereco_step,
    # process.dimurerecowrongsign_step,
    process.p
)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## Add the event selection filters. Going to be used in evtSel ###############################################################
process.Flag_colEvtSel = cms.Path(process.offlinePrimaryVerticesRecovery * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.offlinePrimaryVerticesRecovery * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.offlinePrimaryVerticesRecovery * process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.offlinePrimaryVerticesRecovery * process.clusterCompatibilityFilter)
# eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter , process.Flag_clusterCompatibilityFilter ]
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter ]

#! Adding the process schedule !!!!!!!!!!!!!!!!!!
for P in eventFilterPaths:
    process.schedule.insert(0, P)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!