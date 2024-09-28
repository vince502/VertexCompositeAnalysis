diff --git a/VertexCompositeProducer/test/pPbSkimAndTree2016_DoubleD0Both_BDT_MC_ForLamC_cfg.py b/VertexCompositeProducer/test/pPbSkimAndTree2016_DoubleD0Both_BDT_MC_ForLamC_cfg.py
index 620e283..d9029ea 100644
--- a/VertexCompositeProducer/test/pPbSkimAndTree2016_DoubleD0Both_BDT_MC_ForLamC_cfg.py
+++ b/VertexCompositeProducer/test/pPbSkimAndTree2016_DoubleD0Both_BDT_MC_ForLamC_cfg.py
@@ -32,7 +32,6 @@ process.source = cms.Source("PoolSource",
 #'file:682C6215-7F9A-E711-AB29-0025901D49AC.root'
 #'file:/eos/home-s/soohwan/Analysis/DmesonpPb/MCGen/EvtGenMod/CMSSW_8_0_30/src/HIN-pPb816Summer16DR-00033.root'
 'file:HIN-pPb816Summer16DR-00033_68.root',
-'file:HIN-pPb816Summer16DR-00033_97.root',
 #'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_1.root',
 #'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_10.root',
 #'file:/eos/home-s/soohwan/Analysis/DmesonpPb/VertexCompositeTree/CMSSW_8_0_36_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/MCAOD_DD_Part/HIN-pPb816Summer16DR-00033_11.root',
@@ -137,6 +136,7 @@ process.generalD0CandidatesNewWrongSign = process.generalD0CandidatesNew.clone(i
 
 process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDDCandidates_cff")
 process.generalDDCandidatesNew = process.generalDDCandidates.clone()
+process.generalDDCandidatesNew.d0Collection= cms.InputTag("d0selectorMCNewReduced")
 process.generalDDCandidatesNew.trkPtSumCut = cms.double(0.0)
 process.generalDDCandidatesNew.trkEtaDiffCut = cms.double(999.0)
 process.generalDDCandidatesNew.tkNhitsCut = cms.int32(0)
@@ -190,13 +190,22 @@ process.d0selectorMCNewReduced.trkEtaDiffMax = cms.untracked.double(1.)
 process.d0selectorMCNewReduced.trkNHitMin = cms.untracked.int32(11)
 process.d0selectorMCNewReduced.cand3DPointingAngleMax = cms.untracked.double(1.0)
 process.d0selectorMCNewReduced.cand2DPointingAngleMax = cms.untracked.double(1.0)
-process.d0selectorMCNewReduced.candpTMin = cms.untracked.double(2.0)
-process.d0selectorMCNewReduced.candYMin = cms.untracked.double(-1.0)
-process.d0selectorMCNewReduced.candYMax = cms.untracked.double(1.0)
+process.d0selectorMCNewReduced.candpTMin = cms.untracked.double(1.95)
+process.d0selectorMCNewReduced.candYMin = cms.untracked.double(-1.1)
+process.d0selectorMCNewReduced.candYMax = cms.untracked.double(1.1)
 process.d0selectorMCNewReduced.mvaMin = cms.untracked.double(-1)
 
-process.generalDDCandidatesNew.d0Collection = cms.InputTag("d0selectorMCNewReduced:D0")
-process.generalDDCandidatesNew.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")
+#process.d0selectorMCNewReduced.trkPtMin = cms.untracked.double(0.0)
+#process.d0selectorMCNewReduced.trkPtSumMin = cms.untracked.double(0.0)
+#process.d0selectorMCNewReduced.trkEtaDiffMax = cms.untracked.double(999.0)
+#process.d0selectorMCNewReduced.trkNHitMin = cms.untracked.int32(0)
+#process.d0selectorMCNewReduced.cand3DPointingAngleMax = cms.untracked.double(999.0)
+#process.d0selectorMCNewReduced.cand2DPointingAngleMax = cms.untracked.double(999.0)
+#process.d0selectorMCNewReduced.candpTMin = cms.untracked.double(0.0)
+#process.d0selectorMCNewReduced.candYMin = cms.untracked.double(5.0)
+#process.d0selectorMCNewReduced.candYMax = cms.untracked.double(5.0)
+#process.d0selectorMCNewReduced.mvaMin = cms.untracked.double(-1)
+
 
 process.d0ana_mc_newreduced = process.d0ana_mc.clone()
 process.d0ana_mc_newreduced.saveTree = True
@@ -204,6 +213,7 @@ process.d0ana_mc_newreduced.VertexCompositeCollection = cms.untracked.InputTag("
 process.d0ana_mc_newreduced.MVACollection = cms.InputTag("d0selectorMCNewReduced:MVAValuesNewD0")
 process.d0ana_mc_newreduced.DCAValCollection = cms.InputTag("d0selectorMCNewReduced:DCAValuesNewD0")
 process.d0ana_mc_newreduced.DCAErrCollection = cms.InputTag("d0selectorMCNewReduced:DCAErrorsNewD0")
+process.d0ana_mc_newreduced.MVACollection2= cms.InputTag("")
 
 process.ddana_new = process.ddana_mc.clone()
 process.ddana_new.twoLayerDecay = cms.untracked.bool(True)
@@ -229,6 +239,7 @@ process.ddana_old.MVACollection2= cms.InputTag("generalDDCandidatesNew:MVAValues
 
 
 process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorMCNewReduced * process.d0ana_mc_newreduced * process.generalDDCandidatesNew * process.ddana_new )
+#process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorMCNewReduced)# * process.d0ana_mc_newreduced * process.generalDDCandidatesNew * process.ddana_new )
 
 # eventinfoana must be in EndPath, and process.eventinfoana.selectEvents must be the name of eventFilter_HM Path
 process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
@@ -278,12 +289,12 @@ process.schedule = cms.Schedule(
 # for P in eventFilterPaths:
 #     process.schedule.insert(0, P)
 
-#process.output = cms.OutputModule("PoolOutputModule",
-#    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM"),
-#    fileName = cms.untracked.string('output.root'),
-#)
+process.output = cms.OutputModule("PoolOutputModule",
+    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM"),
+    fileName = cms.untracked.string('output.root'),
+)
 
-#process.outputPath = cms.EndPath(process.output)
-#process.schedule.append(process.outputPath)
+process.outputPath = cms.EndPath(process.output)
+process.schedule.append(process.outputPath)
 
 #process.options.numberOfThreads = cms.untracked.int32(1)
