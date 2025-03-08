import FWCore.ParameterSet.Config as cms

process = cms.Process("writeGBRForests")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1) # NB: needs to be set to 1 so that GBRForestWriter::analyze method gets called exactly once
)

process.source = cms.Source("EmptySource")

process.load('Configuration/StandardSequences/Services_cff')

process.gbrForestWriter = cms.EDAnalyzer("GBRForestWriter",
    jobs = cms.VPSet(
        cms.PSet(
            inputFileName = cms.FileInPath('VertexCompositeAnalysis/VertexCompositeProducer/data/xgboost.xml'),
            inputFileType = cms.string("XML"),
            #inputVariables = cms.vstring( 'VtxProb', 'dca3D', 'v3DCosPointingAngle', 'v3DPointingAngle', 'v2DCosPointingAngle', 'v2DPointingAngle', 'v3DDecayLengthSignificance', 'v3DDecayLength', 'v2DDecayLengthSignificance', 'v2DDecayLength', 'pTD1', 'EtaD1', 'pTD2', 'EtaD2' ),
            inputVariables = cms.vstring( 
		'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11','f12','f13','f14',
	     ),
            spectatorVariables = cms.vstring(),
            methodName = cms.string("BDT"),
            gbrForestName = cms.string("DStarInPbPb"),
            outputFileType = cms.string("GBRForest"),
            outputFileName = cms.string("GBRForestfile_XGBDT_PromptDstarInPbPb_default_MB_OnlyMC_v3.root")
        )
    )
)



process.p = cms.Path(process.gbrForestWriter)
