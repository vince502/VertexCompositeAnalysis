import FWCore.ParameterSet.Config as cms

# D0 -> Kπππ (4-track) analyzer configuration
# Mirrors the structure of the standard D0 analyzer but points to the 4-track producer outputs.
d04pana = cms.EDAnalyzer(
    'PATCompositeTreeProducer3',
    doRecoNtuple=cms.untracked.bool(True),
    doGenNtuple=cms.untracked.bool(False),
    doGenMatching=cms.untracked.bool(False),
    doGenMatchingTOF=cms.untracked.bool(False),
    hasSwap=cms.untracked.bool(False),
    decayInGen=cms.untracked.bool(False),
    twoLayerDecay=cms.untracked.bool(False),
    threeProngDecay=cms.untracked.bool(False),
    PID=cms.untracked.int32(421),
    PID_dau1=cms.untracked.int32(211),
    PID_dau2=cms.untracked.int32(321),
    deltaR=cms.untracked.double(0.03),
    VertexCollection=cms.untracked.InputTag('offlinePrimaryVertices'),
    TrackCollection=cms.untracked.InputTag('generalTracks'),
    CompositeCollection=cms.untracked.InputTag('generalD04PCandidates:D04P'),
    GenParticleCollection=cms.untracked.InputTag('genParticles'),
    MuonCollection=cms.untracked.InputTag('null'),
    doMuon=cms.untracked.bool(False),
    doMuonFull=cms.untracked.bool(False),
    saveTree=cms.untracked.bool(True),
    saveHistogram=cms.untracked.bool(False),
    saveAllHistogram=cms.untracked.bool(False),
    massHistPeak=cms.untracked.double(1.865),
    massHistWidth=cms.untracked.double(0.3),
    massHistBins=cms.untracked.int32(120),
    pTBins=cms.untracked.vdouble(0, 1.2, 1.5, 2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0),
    yBins=cms.untracked.vdouble(-2.4, -1.0, 0.0, 1.0, 2.4),
    useAnyMVA=cms.bool(False),
    isSkimMVA=cms.untracked.bool(False),
    MVACollection=cms.InputTag('generalD04PCandidates:MVAValuesD04P'),
    isCentrality=cms.bool(False),
    centralityBinLabel=cms.InputTag('centralityBin', 'HFtowers'),
    centralitySrc=cms.InputTag('hiCentrality')
)

d04pana_mc = d04pana.clone(
    doGenNtuple=cms.untracked.bool(True),
    doGenMatching=cms.untracked.bool(True),
    decayInGen=cms.untracked.bool(True),
    saveAllHistogram=cms.untracked.bool(True)
)
