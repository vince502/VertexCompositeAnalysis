import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/SKIM/SKIM_AOD_HIForward0_HIRun2023A_20231009/HIForward0/SKIM_AOD_HIForward0_HIRun2023A_20231009/231009_081732/0000/reco_RAW2DIGI_L1Reco_RECO_89.root"),
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/dileptons/Data2023/DimuonSkims/DimuonSkim_PhysicsHIPhysicsRawPrime0_Run374810_ls140to215_miniAOD.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v4')

## Set ZDC information
#process.es_pool = cms.ESSource("PoolDBESSource",
#    timetype = cms.string('runnumber'),
#    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
#    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#    authenticationMethod = cms.untracked.uint32(1)
#)
#process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
#process.es_ascii = cms.ESSource('HcalTextCalibrations',
#    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("emap_2023_newZDC_v3.txt")))
#)

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db"),
        label = cms.untracked.string("HFtowers")
        )
    ]
)
process.cent_seq = cms.Sequence(process.centralityBin)

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

# DiMu selection
SoftIdReco = "(innerTrack.isNonnull && innerTrack.quality(\"highPurity\") && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && isTrackerMuon)"
muonSelection = cms.string("(pt > 0.0 && abs(eta) < 2.5) && "+SoftIdReco)
diMuSelection = cms.string("charge==0")
hpMuSelection = cms.string("innerTrack.isNonnull && innerTrack.quality(\"highPurity\")")
process.diMu = generalParticles.clone(
    pdgId = cms.int32(443),
    preSelection = diMuSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(13), charge = cms.int32(+1), selection = muonSelection),
        cms.PSet(pdgId = cms.int32(13), charge = cms.int32(-1), selection = muonSelection),
    ]),
)
process.oneDiMu = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diMu"), minNumber = cms.uint32(1))

# Add muons
#from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
#doPATMuons(process)


#from RecoMuon.MuonIdentification.calomuons_cfi import calomuons
#process.mergedMuons = cms.EDProducer("CaloMuonMerger",
#    muons    = cms.InputTag("muons::RECO"),
#    caloMuons = cms.InputTag("calomuons"),
#    minCaloCompatibility = calomuons.minCaloCompatibility,
#    mergeTracks = cms.bool(True),
#    tracks = cms.InputTag("generalTracks"),
#    mergeCaloMuons = cms.bool(False),
#    caloMuonsCut = cms.string(""),
#    muonsCut     = hpMuSelection,
#    tracksCut    = cms.string("quality(\"highPurity\")"),
#)

# Add diMu event selection
process.twoMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"), minNumber = cms.uint32(2))
process.hpMuons = cms.EDFilter("MuonSelector", src = cms.InputTag("muons"), cut = hpMuSelection)
process.maxTwoHPMuons = cms.EDFilter("PATCandViewCountFilter", src = cms.InputTag("hpMuons"), minNumber = cms.uint32(0), maxNumber = cms.uint32(2))
process.goodMuons = cms.EDFilter("MuonSelector",
            src = cms.InputTag("muons"),
            cut = muonSelection,
            )
process.twoGoodMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodMuons"), minNumber = cms.uint32(2))
process.goodDiMuons = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = diMuSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodMuons@+ goodMuons@-')
            )
process.oneGoodDiMu = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDiMuons"), minNumber = cms.uint32(1))
process.diMuEvtSel = cms.Sequence( process.twoMuons * process.hpMuons * process.maxTwoHPMuons * process.goodMuons * process.twoGoodMuons * process.goodDiMuons * process.oneGoodDiMu)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    "HLT_HIL*SingleMu*_v*", "HLT_HIL*DoubleMu*_v*", "HLT_HIMinimumBiasHF1AND*_v*"
]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter)
process.primaryVertexFilter.src = cms.InputTag("offlineSlimmedPrimaryVertices")
process.primaryVertexFilter.cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2")

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.diMuEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedTracksAndVertices_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedMuons_cfi')
process.diMu_rereco_step = cms.Path(process.eventFilter_HM  * process.unpackedTracksAndVertices * process.unpackedMuons  *process.patMuonSequence *  process.diMu * process.oneDiMu * process.cent_seq)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.diMuAna = particleAna.clone(
  recoParticles = cms.InputTag("diMu"),
  selectEvents = cms.string("diMu_rereco_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
  ),
  triggerInfo = cms.untracked.VPSet([
    # UPC muon triggers
    cms.PSet(path = cms.string('HLT_HIL*SingleMu*_v*')),
    cms.PSet(path = cms.string('HLT_HIL*DoubleMu*_v*')),
    cms.PSet(path = cms.string('HLT_HIMinimumBiasHF1AND*_v*')),
  ]),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('diMu_ana.root'))
process.p = cms.EndPath(process.diMuAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.diMu_rereco_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)

eventFilterPaths = [ 
    process.Flag_colEvtSel ,
    process.Flag_hfCoincFilter2Th4 ,
    process.Flag_primaryVertexFilter,
]

process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.primaryVertexFilter *
    process.diMuEvtSel
)

for P in eventFilterPaths:
    process.schedule.insert(0, P)

# Add recovery for offline primary vertex
from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
process = MassReplaceInputTag(process, "muons", "slimmedMuons")
process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")
process = MassReplaceInputTag(process,"genParticles","prunedGenParticles")
#process.mergedMuons.muons = cms.InputTag("muons")
