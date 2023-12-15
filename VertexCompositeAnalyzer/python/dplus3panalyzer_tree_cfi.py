import FWCore.ParameterSet.Config as cms

dplus3pana = cms.EDAnalyzer('VertexCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(411),
  PID_dau1 = cms.untracked.int32(421),
  PID_dau2 = cms.untracked.int32(211),
  PID_dau3 = cms.untracked.int32(211),
  deltaR = cms.untracked.double(0.03),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDPlus3PCandidatesNew:DPlus3P"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),
  doMuonFull = cms.untracked.bool(False),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),
  massHistPeak = cms.untracked.double(2.01),
  massHistWidth = cms.untracked.double(0.2),
  massHistBins = cms.untracked.int32(100),

  pTBins = cms.untracked.vdouble(0,1.2,1.5,2.4,3.0,3.5,4.2,5.0,6.0,7.0,8.0),
  yBins = cms.untracked.vdouble(-2.4,-1.0,0.0,1.0,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  #MVACollection = cms.InputTag("generalD0CandidatesNew:MVAValues"),

  isCentrality = cms.bool(False),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality")
                              )

dplus3pana_mc = cms.EDAnalyzer('VertexCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(411),
  PID_dau1 = cms.untracked.int32(421),
  PID_dau2 = cms.untracked.int32(211),
  PID_dau3 = cms.untracked.int32(211),
  deltaR = cms.untracked.double(0.03),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDPlus3PCandidatesNew:DPlus3P"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),
  doMuonFull = cms.untracked.bool(False),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(True),
  massHistPeak = cms.untracked.double(1.86),
  massHistWidth = cms.untracked.double(0.2),
  massHistBins = cms.untracked.int32(100),

  pTBins = cms.untracked.vdouble(0,1.2,1.5,2.4,3.0,3.5,4.2,5.0,6.0,7.0,8.0),
  yBins = cms.untracked.vdouble(-2.4,-1.0,0.0,1.0,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalDPlus3PCandidatesNew:MVAValues")
                              )
