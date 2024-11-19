import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras

#----------------------------------------------------------------------------

# Setup Settings for ONIA TREE: 2024 pp ref data

HLTProcess     = "HLT" # Name of HLT process
isMC           = False # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "GlbTrk" # Single muon selection: All, Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker), GlbOrTrk, TwoGlbAmongThree (which requires two isGlobal for a trimuon, and one isGlobal for a dimuon) are available
applyEventSel  = True # Only apply Event Selection if the required collections are present
OnlySoftMuons  = False # Keep only isSoftMuon's (without highPurity, and without isGlobal which should be put in 'muonSelection' parameter) from the beginning of HiSkim. If you want the full SoftMuon selection, set this flag false and add 'isSoftMuon' in lowerPuritySelection. In any case, if applyCuts=True, isSoftMuon is required at HiAnalysis level for muons of selected dimuons.
applyCuts      = False # At HiAnalysis level, apply kinematic acceptance cuts + identification cuts (isSoftMuon (without highPurity) or isTightMuon, depending on TightGlobalMuon flag) for muons from selected di(tri)muons + hard-coded cuts on the di(tri)muon that you would want to add (but recommended to add everything in LateDimuonSelection, applied at the end of HiSkim)
SumETvariables = False  # Whether to write out SumET-related variables
SofterSgMuAcceptance = True # Whether to accept muons with a softer acceptance cuts than the usual (pt>3.5GeV at central eta, pt>1.5 at high |eta|). Applies when applyCuts=True
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
print( "[INFO] OneMatchedHLTMu      = " + (OneMatchedHLTMu if OneMatchedHLTMu > -1 else "False") )
print( "[INFO] miniAOD              = " + ("True" if miniAOD else "False") )
print( "[INFO] UsePropToMuonSt      = " + ("True" if UsePropToMuonSt else "False") )
print( " " )

# set up process
process = cms.Process("HIOnia", eras.Run3_2024_ppRef)

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Name
options.outputFile = "Oniatree_ppData2024_miniAOD.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles =[
#  '/store/data/Run2023F/PPRefDoubleMuon1/MINIAOD/PromptReco-v1/000/373/710/00000/1f07f64a-12c6-41a8-9440-317c201cb3b1.root',
#  '/store/data/Run2023F/PPRefDoubleMuon1/MINIAOD/PromptReco-v1/000/373/710/00000/b13150ae-54e0-4a80-8d7f-2203a6773140.root',
#  '/store/data/Run2023F/PPRefDoubleMuon3/MINIAOD/PromptReco-v1/000/373/710/00000/b72132ca-d17b-4682-a211-996dfe06b528.root',
#  '/store/data/Run2023F/PPRefDoubleMuon3/MINIAOD/PromptReco-v1/000/373/710/00000/d5d569e1-632d-4991-87c4-d67571def74e.root',
#  '/store/data/Run2023F/PPRefDoubleMuon2/MINIAOD/PromptReco-v1/000/373/710/00000/95ac62f6-da92-4f81-ad5b-7aaeb6dcca0c.root'
#'file:/eos/cms/store/express/Run2024J/ExpressPhysics/FEVT/Express-v1/000/387/474/00000/2f267e89-8367-44b2-8442-58e4dfd05aa1.root',
#'file:/afs/cern.ch/work/s/soohwan/private/Analysis/General2024Analysis/CMSSW_14_1_4_patch1/src/cfg_txt/recoppraw2mini_RAW2DIGI_L1Reco_RECO_PAT.root',
#'file:/afs/cern.ch/work/s/soohwan/private/Analysis/General2024Analysis/CMSSW_14_1_4_patch1/src/cfg_txt/recoppraw2mini_RAW2DIGI_L1Reco_RECO_PAT.root',
#'/store/data/Run2024J/PPRefDoubleMuon0/MINIAOD/PromptReco-v1/000/387/695/00000/3be8f20e-6df1-4678-baed-b3c65b5ac756.root',
'file:3be8f20e-6df1-4678-baed-b3c65b5ac756.root',
#'/store/hidata/HIRun2023A/HIForward1/MINIAOD/16Jan2024-v1/40000/574b9ae2-e2d7-45db-96af-ab91bfac097e.root',

]
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

triggerList    = {
		# Double Muon Trigger List
		'DoubleMuonTrigger' : cms.vstring(
			      "HLT_PPRefL1DoubleMu0_Open_v",
            "HLT_PPRefL1DoubleMu0_v",
            "HLT_PPRefL1DoubleMu0_SQ_v",
            "HLT_PPRefL1DoubleMu2_v",
            "HLT_PPRefL1DoubleMu2_SQ_v",
            "HLT_PPRefL2DoubleMu0_Open_v",
            "HLT_PPRefL2DoubleMu0_v",
            "HLT_PPRefL3DoubleMu0_Open_v",
            "HLT_PPRefL3DoubleMu0_v"
            ),
        # Single Muon Trigger List
        'SingleMuonTrigger' : cms.vstring(
            "HLT_PPRefL1SingleMu7_v",
            "HLT_PPRefL1SingleMu12_v",
            "HLT_PPRefL2SingleMu7_v",
            "HLT_PPRefL2SingleMu12_v",
            "HLT_PPRefL2SingleMu15_v",
            "HLT_PPRefL2SingleMu20_v",
            "HLT_PPRefL3SingleMu3_v",
            "HLT_PPRefL3SingleMu5_v",
            "HLT_PPRefL3SingleMu7_v",
            "HLT_PPRefL3SingleMu12_v",
            "HLT_PPRefL3SingleMu15_v",
            "HLT_PPRefL3SingleMu20_v",
			)
                }

# Global tag, see https://github.com/cms-sw/cmssw/blob/master/Configuration/AlCa/python/autoCond.py
if isMC:
  globalTag = 'auto:phase1_2024_realistic_ppRef5TeV' #for Run3 MC : phase1_2023_realistic
else:
  globalTag = '141X_dataRun3_Prompt_v3'

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


#----------------------------------------------------------------------------

# For OniaTree Analyzer
from HiAnalysis.HiOnia.oniaTreeAnalyzer_cff import oniaTreeAnalyzer
oniaTreeAnalyzer(process,
                 muonTriggerList=triggerList, #HLTProName=HLTProcess,
                 muonSelection=muonSelection, L1Stage=2, isMC=isMC, pdgID=pdgId, outputFileName=options.outputFile, doTrimu=doTrimuons#, OnlySingleMuons=True
)

process.onia2MuMuPatGlbGlb.dimuonSelection       = cms.string("mass > 2 && charge==0 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
#process.onia2MuMuPatGlbGlb.lowerPuritySelection  = cms.string("pt > 5 || isPFMuon || (pt>1.2 && (isGlobalMuon || isStandAloneMuon)) || (isTrackerMuon && track.quality('highPurity'))")
#process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("") ## No need to repeat lowerPuritySelection in there, already included
if applyCuts:
  process.onia2MuMuPatGlbGlb.LateDimuonSel         = cms.string("userFloat(\"vProb\")>0.01")
process.onia2MuMuPatGlbGlb.onlySoftMuons         = cms.bool(OnlySoftMuons)
process.hionia.minimumFlag      = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.useGeTracks      = cms.untracked.bool(keepExtraColl) #for Reco_trk_*
process.hionia.fillRecoTracks   = cms.bool(keepExtraColl)           #for Reco_trk_*
#process.hionia.CentralitySrc    = cms.InputTag("hiCentrality")
#process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
#process.hionia.muonLessPV       = cms.bool(False)
process.hionia.SofterSgMuAcceptance = cms.bool(SofterSgMuAcceptance)
process.hionia.SumETvariables   = cms.bool(SumETvariables)
process.hionia.applyCuts        = cms.bool(applyCuts)
process.hionia.AtLeastOneCand   = cms.bool(atLeastOneCand)
process.hionia.OneMatchedHLTMu  = cms.int32(OneMatchedHLTMu)
process.hionia.checkTrigNames   = cms.bool(False)#change this to get the event-level trigger info in hStats output (but creates lots of warnings when fake trigger names are used)
process.hionia.mom4format       = cms.string(useMomFormat)
process.hionia.isHI = cms.untracked.bool(False)
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



if applyEventSel:
    # Offline event filters
    process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')

    # HLT trigger firing events
    import HLTrigger.HLTfilters.hltHighLevel_cfi
    process.hltHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
    process.hltHI.HLTPaths = ["HLT_PPRefL*DoubleMu*_v*","HLT_PPRefL*SingleMu*_v*","HLT_HIUPC*Mu*v*"]
    process.hltHI.throw = False
    process.hltHI.andOr = True
    process.oniaTreeAna.replace(process.patMuonSequence, process.primaryVertexFilter * process.hltHI * process.patMuonSequence )

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

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalOttCandidates_cff")
process.generalOttCandidatesNew = process.generalOttCandidates.clone()
process.generalOttCandidatesNew.dimuons = cms.InputTag('onia2MuMuPatGlbGlb')
process.generalOttCandidatesNew.vertexRecoAlgorithm = cms.InputTag('unpackedTracksAndVertices')
process.generalOttCandidatesNew.trackRecoAlgorithm = cms.InputTag('unpackedTracksAndVertices')

#process.generalOttCandidatesNew.vertexRecoAlgorithm = cms.InputTag('offlineSlimmedPrimaryVertices')
process.generalOttCandidatesNew.usePixelTracks = cms.bool(False)
process.generalOttCandidatesNew.pixelTracks = cms.InputTag('unpackedPixelTracks')


process.generalOttCandidatesNew.batTrkPtSumCut = cms.double(0.1)
process.generalOttCandidatesNew.batTrkEtaDiffCut = cms.double(1.6)
process.generalOttCandidatesNew.batTkChi2Cut = cms.double(3)
process.generalOttCandidatesNew.batTkNhitsCut = cms.int32(0)
process.generalOttCandidatesNew.batTkPtErrCut = cms.double(0.1)
process.generalOttCandidatesNew.batTkPtCut = cms.double(0.10)
process.generalOttCandidatesNew.alphaCut = cms.double(999.0)
process.generalOttCandidatesNew.alpha2DCut = cms.double(999.0)
process.generalOttCandidatesNew.bPtCut = cms.double(0.0)
process.generalOttCandidatesNew.bVtxChiProbCut = cms.double(0.010)
process.generalOttCandidatesNew.mPiKCutMin = cms.double(0.0)
process.generalOttCandidatesNew.mPiKCutMax = cms.double(40.0)
process.generalOttCandidatesNew.bMassCut = cms.double(12)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.ottanalyzer_tree_cff")
process.ottana_new = process.ottana.clone()
process.ottana_new.VertexCollection = cms.untracked.InputTag('unpackedTracksAndVertices')
process.ottana_new.TrackCollection = cms.untracked.InputTag('unpackedTracksAndVertices')
process.ottana_new.doRecoNtuple = True
process.ottana_new.PID = 20443

process.ottana_new.threeProngDecay = False 
process.ottana_new.PID_dau1 = 443
process.ottana_new.PID_dau2 = 113

#process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
#process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 140X, mc")
process.load("HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_ppref_data_cff")
# use data version to avoid PbPb MC
#process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
#process.hiEvtAnalyzer.Vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.hiEvtAnalyzer.doCentrality = cms.bool(False)
#process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
#process.hiEvtAnalyzer.doEvtPlaneFlat = cms.bool(False)
#process.hiEvtAnalyzer.doMC = cms.bool(True) # general MC info
#process.hiEvtAnalyzer.doHiMC = cms.bool(False) # HI specific MC info
#process.hiEvtAnalyzer.doHFfilters = cms.bool(False) # Disable HF filters for ppRef


#process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
#process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
#process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

## TODO: Many of these triggers are not available in the test file
#from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
#process.hltobject.triggerNames = trigger_list_mc

# Gen particles
#process.load('HeavyIonsAnalysis.EventAnalysis.HiGenAnalyzer_cfi')

process.forest = cms.Sequence(
    process.primaryVertexFilter
#    process.HiForestInfo
#    process.hltanalysis *
#    process.hiEvtAnalyzer *
#    process.hltobject +
#    process.l1object +
#    process.HiGenParticleAna
)

addR3Jets = False
addR4Jets = True
addR5Jets = True

if addR3Jets or addR4Jets :
#    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupPprefJets

    if addR3Jets :
        process.jetsR3 = cms.Sequence()
        setupPprefJets('ak3PF', process.jetsR3, process, isMC = 0, radius = 0.30, JECTag = 'AK3PF')
        process.ak3PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.ak3PFJetAnalyzer = process.ak4PFJetAnalyzer.clone(jetTag = "ak3PFpatJets", jetName = 'ak3PF', genjetTag = "ak3GenJetsNoNu")
        process.forest += process.jetsR3 * process.ak3PFJetAnalyzer

    if addR4Jets :
        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections
        process.jetsR4 = cms.Sequence()
        setupPprefJets('ak04PF', process.jetsR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF')
        process.ak04PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.ak04PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.ak4PFJetAnalyzer.jetTag = 'ak04PFpatJets'
        process.ak4PFJetAnalyzer.jetName = 'ak04PF'
        process.ak4PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents.
        process.forest += process.jetsR4 * process.ak4PFJetAnalyzer
#    if addR5Jets :
#        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections
#        process.jetsR5 = cms.Sequence()
#        setupPprefJets('ak05PF', process.jetsR5, process, isMC = 0, radius = 0.50, JECTag = 'AK5PF')
#        process.ak05PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
#        process.ak05PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
#        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
#        process.ak5PFJetAnalyzer.jetTag = 'ak05PFpatJets'
#        process.ak5PFJetAnalyzer.jetName = 'ak05PF'
#        process.ak5PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents.
#        process.forest += process.jetsR5 * process.ak5PFJetAnalyzer


#process.schedule.append( process.forest )

#process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter  * process.generalOttCandidatesNew * process.forest * process.ottana_new)
process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter  * process.generalOttCandidatesNew * process.ottana_new)
#process.generalOttCandidatesNew.d0RecoAlgorithm = cms.InputTag('')

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

process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.schedule  = cms.Schedule( process.oniaTreeAna )

#process.output = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring(["drop *", "keep *_onia2MuMu*_*_HIOnia", "keep *_generalOtt*_*_HIOnia"]),
#    #fileName = cms.untracked.string("/eos/cms/store/group/phys_heavyions/soohwan/Run3_2024/ppRef_2024/output2.root"),
#    fileName = cms.untracked.string("output2.root"),
#)
#process.output_path = cms.EndPath(process.output)
#process.schedule.append( process.output_path )
