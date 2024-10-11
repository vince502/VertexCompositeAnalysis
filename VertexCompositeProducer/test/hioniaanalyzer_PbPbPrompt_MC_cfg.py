import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras

#----------------------------------------------------------------------------

# Setup Settings for ONIA TREE: 2024 PbPb MC

HLTProcess     = "HLT" # Name of HLT process
isMC           = True # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "All" # Single muon selection: All, Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker), GlbOrTrk, TwoGlbAmongThree (which requires two isGlobal for a trimuon, and one isGlobal for a dimuon) are available
applyEventSel  = True # Only apply Event Selection if the required collections are present
OnlySoftMuons  = False # Keep only isSoftMuon's (without highPurity, and without isGlobal which should be put in 'muonSelection' parameter) from the beginning of HiSkim. If you want the full SoftMuon selection, set this flag false and add 'isSoftMuon' in lowerPuritySelection. In any case, if applyCuts=True, isSoftMuon is required at HiAnalysis level for muons of selected dimuons.
applyCuts      = False # At HiAnalysis level, apply kinematic acceptance cuts + identification cuts (isSoftMuon (without highPurity) or isTightMuon, depending on TightGlobalMuon flag) for muons from selected di(tri)muons + hard-coded cuts on the di(tri)muon that you would want to add (but recommended to add everything in LateDimuonSelection, applied at the end of HiSkim)
SumETvariables = True  # Whether to write out SumET-related variables
SofterSgMuAcceptance = False # Whether to accept muons with a softer acceptance cuts than the usual (pt>3.5GeV at central eta, pt>1.5 at high |eta|). Applies when applyCuts=True
doTrimuons     = False # Make collections of trimuon candidates in addition to dimuons, and keep only events with >0 trimuons (if atLeastOneCand)
doDimuonTrk    = False # Make collections of Jpsi+track candidates in addition to dimuons
atLeastOneCand = False # Keep only events that have one selected dimuon (or at least one trimuon if doTrimuons = true). BEWARE this can cause trouble in .root output if no event is selected by onia2MuMuPatGlbGlbFilter!
OneMatchedHLTMu = -1   # Keep only di(tri)muons of which the one(two) muon(s) are matched to the HLT Filter of this number. You can get the desired number in the output of oniaTree. Set to -1 for no matching.
#############################################################################
keepExtraColl  = False # General Tracks + Stand Alone Muons + Converted Photon collections
miniAOD        = True # whether the input file is in miniAOD format (default is AOD)
UsePropToMuonSt = True # whether to use L1 propagated muons (works only for miniAOD now)
pdgId = 443 # J/Psi : 443, Y(1S) : 553
useMomFormat = "vector" # default "array" for TClonesArray of TLorentzVector. Use "vector" for std::vector<float> of pt, eta, phi, M
#----------------------------------------------------------------------------

# Print Onia Tree settings:
print( " " )
print( "[INFO] Settings used for ONIA TREE: " )
print( "[INFO] isMC                 = " + ("True" if isMC else "False") )
print( "[INFO] applyEventSel        = " + ("True" if applyEventSel else "False") )
print( "[INFO] applyCuts            = " + ("True" if applyCuts else "False") )
print( "[INFO] keepExtraColl        = " + ("True" if keepExtraColl else "False") )
print( "[INFO] SumETvariables       = " + ("True" if SumETvariables else "False") )
print( "[INFO] SofterSgMuAcceptance = " + ("True" if SofterSgMuAcceptance else "False") )
print( "[INFO] muonSelection        = " + muonSelection )
print( "[INFO] onlySoftMuons        = " + ("True" if OnlySoftMuons else "False") )
print( "[INFO] doTrimuons           = " + ("True" if doTrimuons else "False") )
print( "[INFO] doDimuonTrk          = " + ("True" if doDimuonTrk else "False") )
print( "[INFO] atLeastOneCand       = " + ("True" if atLeastOneCand else "False") )
print( "[INFO] OneMatchedHLTMu      = " + ("True" if OneMatchedHLTMu > -1 else "False") )
print( "[INFO] miniAOD              = " + ("True" if miniAOD else "False") )
print( "[INFO] UsePropToMuonSt      = " + ("True" if UsePropToMuonSt else "False") )
print( " " )

# set up process
process = cms.Process("HIOnia", eras.Run3_pp_on_PbPb)

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Name
options.outputFile = "Oniatree_MC_miniAOD.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles =[
  # '/store/mc/HINPbPbSpring23MiniAOD/PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/MINIAODSIM/132X_mcRun3_2023_realistic_HI_v9-v2/130000/0c311838-da71-4010-b4f9-313c7fef7745.root',
  "file:/eos/cms/store/group/phys_heavyions/dileptons/Data2023/MINIAOD/HIPhysicsRawPrime0/Run375064/7ed5766f-6b1d-415e-8916-e62825a6347f.root",
]
options.maxEvents = 100 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

triggerList    = {
		# Double Muon Trigger List
		'DoubleMuonTrigger' : cms.vstring(
                        "HLT_HIL1DoubleMu0_MaxDr3p5_Open_v",#0
                        "HLT_HIL1DoubleMu0_v",#1
                        "HLT_HIL1DoubleMu0_SQ_v",#2
                        "HLT_HIL2DoubleMu0_Open_v",#3
                        "HLT_HIL2DoubleMu0_M1p5to6_Open_v",#4
                        "HLT_HIL2DoubleMu2p8_M1p5to6_Open_v",#5
                        "HLT_HIL2DoubleMu0_M7to15_Open_v",#6
                        "HLT_HIL2DoubleMu3_M7to15_Open_v",#7
                        "HLT_HIL3DoubleMu0_M0toInf_Open_v",#8
                        "HLT_HIL3DoubleMu0_Quarkonia_Open_v",#9
                        "HLT_HIL3DoubleMu2_Quarkonia_Open_v",#10
                        "HLT_HIL3DoubleMu0_M2to4p5_Open_v",#11
                        "HLT_HIL3DoubleMu2_M2to4p5_Open_v",#12
                        "HLT_HIL3DoubleMu0_M7to15_Open_v",#13
                        "HLT_HIL3DoubleMu2_M7to15_Open_v",#14
                        ),
                # Single Muon Trigger List
                'SingleMuonTrigger' : cms.vstring(
                        "HLT_HIL1SingleMu0_Open_v",#15
                        "HLT_HIL1SingleMu0_v",#16
                        "HLT_HIL2SingleMu3_Open_v",#17
                        "HLT_HIL2SingleMu5_v",#18
                        "HLT_HIL2SingleMu7_v",#19
                        "HLT_HIL3SingleMu3_Open_v",#20
                        "HLT_HIL3SingleMu5_v",#21
                        "HLT_HIL3SingleMu7_v",#22
                        "HLT_HIL3SingleMu12_v",#23
                        "HLT_HIMinimumBiasHF1AND_v", #24
                        "HLT_HIMinimumBiasHF1ANDZDC2nOR_v", #25
                        "HLT_HIMinimumBiasHF1ANDZDC1nOR_v", #26
			)
                }

## Global tag
if isMC:
  globalTag = 'auto:phase1_2023_realistic_hi' #for Run3 MC : phase1_2023_realistic_hi
else:
  globalTag = '132X_dataRun3_Prompt_v7' # 'auto:run3_data_prompt'

#----------------------------------------------------------------------------

# load the Geometry and Magnetic Field
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')
### For Centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
print('\n\033[31m~*~ USING OFFICIAL MC CENTRALITY TABLE FOR PbPb 2023 ~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_Run3v1302x04_Official_MC"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])

#----------------------------------------------------------------------------

# For OniaTree Analyzer
from HiAnalysis.HiOnia.oniaTreeAnalyzer_cff import oniaTreeAnalyzer
oniaTreeAnalyzer(process,
                 muonTriggerList=triggerList, HLTProName=HLTProcess,
                 muonSelection=muonSelection, L1Stage=2, isMC=isMC, pdgID=pdgId, outputFileName=options.outputFile, doTrimu=doTrimuons,
                 OnlySingleMuons=False
)

#process.onia2MuMuPatGlbGlb.dimuonSelection       = cms.string("8 < mass && mass < 14 && charge==0 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
#process.onia2MuMuPatGlbGlb.lowerPuritySelection  = cms.string("pt > 5 || isPFMuon || (pt>1.2 && (isGlobalMuon || isStandAloneMuon)) || (isTrackerMuon && track.quality('highPurity'))")
#process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("") ## No need to repeat lowerPuritySelection in there, already included
if applyCuts:
  process.onia2MuMuPatGlbGlb.LateDimuonSel         = cms.string("userFloat(\"vProb\")>0.01")
process.onia2MuMuPatGlbGlb.onlySoftMuons         = cms.bool(OnlySoftMuons)
process.hionia.minimumFlag      = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.useGeTracks      = cms.untracked.bool(keepExtraColl) #for Reco_trk_*
process.hionia.fillRecoTracks   = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.CentralitySrc    = cms.InputTag("hiCentrality")
process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
#process.hionia.muonLessPV       = cms.bool(False)
process.hionia.SofterSgMuAcceptance = cms.bool(SofterSgMuAcceptance)
process.hionia.SumETvariables   = cms.bool(SumETvariables)
process.hionia.applyCuts        = cms.bool(applyCuts)
process.hionia.AtLeastOneCand   = cms.bool(atLeastOneCand)
process.hionia.OneMatchedHLTMu  = cms.int32(OneMatchedHLTMu)
process.hionia.checkTrigNames   = cms.bool(False)#change this to get the event-level trigger info in hStats output (but creates lots of warnings when fake trigger names are used)
process.hionia.mom4format       = cms.string(useMomFormat)
process.hionia.genealogyInfo    = cms.bool(True)

'''
#----------------------------------------------------------------------------

# For HLTBitAnalyzer
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.HLTProcessName              = HLTProcess
process.hltbitanalysis.hltresults                  = cms.InputTag("TriggerResults","",HLTProcess)
process.hltbitanalysis.l1tAlgBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.l1tExtBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.gObjectMapRecord            = cms.InputTag("hltGtStage2ObjectMap","",HLTProcess)
process.hltbitanalysis.gmtStage2Digis              = cms.string("hltGtStage2Digis")
process.hltbitanalysis.caloStage2Digis             = cms.string("hltGtStage2Digis")
process.hltbitanalysis.UseL1Stage2                 = cms.untracked.bool(True)
process.hltbitanalysis.getPrescales                = cms.untracked.bool(False)
process.hltbitanalysis.getL1InfoFromEventSetup     = cms.untracked.bool(False)
process.hltbitanalysis.UseTFileService             = cms.untracked.bool(True)
process.hltbitanalysis.RunParameters.HistogramFile = cms.untracked.string(options.outputFile)
process.hltbitanalysis.RunParameters.isData        = cms.untracked.bool(not isMC)
process.hltbitanalysis.RunParameters.Monte         = cms.bool(isMC)
process.hltbitanalysis.RunParameters.GenTracks     = cms.bool(False)
if (HLTProcess == "HLT") :
	process.hltbitanalysis.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.gmtStage2Digis    = cms.string("gtStage2Digis")
	process.hltbitanalysis.caloStage2Digis   = cms.string("gtStage2Digis")

##----------------------------------------------------------------------------

# For HLTObject Analyzer
process.load("HeavyIonsAnalysis.EventAnalysis.hltobject_cfi")
process.hltobject.processName = cms.string(HLTProcess)
process.hltobject.treeName = cms.string(options.outputFile)
process.hltobject.loadTriggersFromHLT = cms.untracked.bool(False)
process.hltobject.triggerNames = triggerList['DoubleMuonTrigger'] + triggerList['SingleMuonTrigger']
process.hltobject.triggerResults = cms.InputTag("TriggerResults","",HLTProcess)
process.hltobject.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","",HLTProcess)

if saveHLT:
  process.oniaTreeAna = cms.Path(process.hltbitanalysis * process.hltobject * process.oniaTreeAna )
'''

process.oniaTreeAna.replace(process.hionia, process.centralityBin * process.hionia )

if applyEventSel:
  process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
  process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
  process.oniaTreeAna.replace(process.patMuonSequence, process.phfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter * process.patMuonSequence )

if atLeastOneCand:
  if doTrimuons:
      process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilterTrimu)
      process.oniaTreeAna.replace(process.patMuonSequence, process.filter3mu * process.pseudoDimuonFilterSequence * process.patMuonSequence)
  elif doDimuonTrk:
      process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilterDimutrk)
      process.oniaTreeAna.replace(process.patMuonSequence, process.pseudoDimuonFilterSequence * process.patMuonSequence)
  else:
      process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter)
      #BEWARE, pseudoDimuonFilterSequence asks for opposite-sign dimuon in given mass range. But saves a lot of time by filtering before running PAT muons
      process.oniaTreeAna.replace(process.patMuonSequence, process.pseudoDimuonFilterSequence * process.patMuonSequence)

process.oniaTreeAna = cms.Path(process.oniaTreeAna)
if miniAOD:
  from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import changeToMiniAOD
  changeToMiniAOD(process)
  process.unpackedMuons.addPropToMuonSt = cms.bool(UsePropToMuonSt)

#----------------------------------------------------------------------------
#Options:
process.source = cms.Source("PoolSource",
#process.source = cms.Source("NewEventStreamFileReader", # for streamer data
		fileNames = cms.untracked.vstring( options.inputFiles ),
		)
process.TFileService = cms.Service("TFileService",
		fileName = cms.string( options.outputFile )
		)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.options.numberOfThreads = 4


########## OTT candidate rereco ###############################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalOttCandidates_cff")
process.generalOttCandidatesNew = process.generalOttCandidates.clone()
process.generalOttCandidatesNew.dimuons = cms.InputTag('oniaGlbGlb')
process.generalOttCandidatesNew.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.generalOttCandidatesNew.trackRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')

process.generalOttCandidatesNew.trkPtSumCut = cms.double(0.0)
process.generalOttCandidatesNew.trkEtaDiffCut = cms.double(4.8)
process.generalOttCandidatesNew.tkNhitsCut = cms.int32(0)
process.generalOttCandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalOttCandidatesNew.tkPtCut = cms.double(1.0)
process.generalOttCandidatesNew.alphaCut = cms.double(1.20)
process.generalOttCandidatesNew.alpha2DCut = cms.double(1.20)
process.generalOttCandidatesNew.dPtCut = cms.double(0.0)
process.generalOttCandidatesNew.VtxChiProbCut = cms.double(0.010)
process.generalOttCandidatesNew.mPiKCutMin = cms.double(3.0)
process.generalOttCandidatesNew.mPiKCutMax = cms.double(14.0)


process.Ottreco_step = cms.Path(  process.generalOttCandidatesNew)

process.schedule  = cms.Schedule( process.oniaTreeAna )

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(["drop *", "keep *_*_*_HIOnia"]),
    fileName = cms.untracked.string("output.root"),
)
process.output_path = cms.EndPath(process.output)
process.schedule.append( process.output_path )
