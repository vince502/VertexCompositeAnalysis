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
            inputFileName = cms.FileInPath('VertexCompositeAnalysis/VertexCompositeProducer/data/xgboost_v2.xml'),
            inputFileType = cms.string("XML"),
            #inputVariables = cms.vstring( 'VtxProb', 'dca3D', 'v3DCosPointingAngle', 'v3DPointingAngle', 'v2DCosPointingAngle', 'v2DPointingAngle', 'v3DDecayLengthSignificance', 'v3DDecayLength', 'v2DDecayLengthSignificance', 'v2DDecayLength', 'pTD1', 'EtaD1', 'pTD2', 'EtaD2' ),
            inputVariables = cms.vstring( 
		'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13',
	     ),
            spectatorVariables = cms.vstring(),
            methodName = cms.string("BDT"),
            gbrForestName = cms.string("D0InpPbXGB"),
            outputFileType = cms.string("GBRForest"),
            outputFileName = cms.string("GBRForestfile_XGBDT_PromptD0InpPb_default_MB_OnlyMC_v2.root")
        )
    )
)

#process.load("CondCore.DBCommon.CondDBCommon_cfi")


#process.PoolDBOutputService = cms.Service("PoolDBOutputService",
#    process.CondDBCommon,
#    timetype = cms.untracked.string('runnumber'),
#    toPut = cms.VPSet(
#        cms.PSet(
#            record = cms.string('btag_CombinedMVAv2_BDT_TMVAv420_74X_v1'),
#            tag = cms.string('btag_CombinedMVAv2_BDT_TMVAv420_74X_v1'),
#            label = cms.untracked.string('btag_CombinedMVAv2_BDT')
#        )
#    )
#)


process.p = cms.Path(process.gbrForestWriter)
