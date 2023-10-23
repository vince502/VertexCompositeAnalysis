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
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/STARLIGHT/2023_10_23/STARLIGHT_5p36TeV_2023Run3/QED_elec_STARLIGHT_5p36TeV_2023Run3_HIFwd_RECO_AOD_2023_10_23/231022_221025/0000/STARLIGHT_QED_elec_RECO_437.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_mcRun3_2023_realistic_HI_v5')

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

# DiMu selection
process.diEle = generalParticles.clone(
    pdgId = cms.int32(443),
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(11), charge = cms.int32(+1)),
        cms.PSet(pdgId = cms.int32(11), charge = cms.int32(-1)),
    ]),
)
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATElectrons
doPATElectrons(process)
process.patElectrons.electronSource = cms.InputTag("gedGsfElectrons")

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hiClusterCompatibility * process.primaryVertexFilter)

# Define the analysis steps
process.diEle_rereco_step = cms.Path(process.patElectrons * process.diEle)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
process.diEleAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("diEle"),
  selectEvents = cms.string(""),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_primaryVertexFilter',
  ),
  triggerInfo = cms.untracked.VPSet([
  ]),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('diEle_ana_mc.root'))
process.p = cms.EndPath(process.diEleAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.diEle_rereco_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)

