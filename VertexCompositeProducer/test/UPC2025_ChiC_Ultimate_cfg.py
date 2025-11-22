import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('ANASKIM', eras.Run3_2025)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(1),
)
process.FastTimerService = cms.Service(
    'FastTimerService',
    printEventSummary=cms.untracked.bool(True),
    printRunSummary=cms.untracked.bool(True),
    printJobSummary=cms.untracked.bool(True),
    enableDQM=cms.untracked.bool(False),
)

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
#        '/store/hidata/HIRun2025A/HIForward0/MINIAOD/PromptReco-v1/000/399/540/00000/491ce449-be44-4fe8-a337-f85c90499ea9.root',
'/store/hidata/HIRun2025A/HIForward0/MINIAOD/PromptReco-v1/000/399/655/00000/04672a9e-7153-4e94-8603-98e3eb08fa58.root'
#'/store/hidata/HIRun2024B/HIForward0/MINIAOD/PromptReco-v1/000/388/317/00000/6bc0d1dc-5ec4-41e3-90f1-34a919e98767.root',
    ),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('151X_dataRun3_Prompt_v1')

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_HIUPC_ZDC1n*',
    'HLT_HIUPC_ZeroBias_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
]

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')
process.colEvtSel = cms.Sequence()

# Track count filter with quality cuts: min 2 tracks, max 5 tracks
# Quality cuts: chi2, nhits > 4, npix > 2
from VertexCompositeAnalysis.VertexCompositeProducer.trackCountFilter_cfi import trackCountFilter
process.trackCountFilter = trackCountFilter.clone()
process.trackCountFilter.trackCollection = cms.InputTag('generalTracks')
process.trackCountFilter.minTrackPt = cms.double(0.1)
process.trackCountFilter.maxTrackEta = cms.double(2.4)
process.trackCountFilter.maxTrackNormalizedChi2 = cms.double(10.0)
process.trackCountFilter.minTrackNHits = cms.int32(5)  # nhits > 4 means >= 5
process.trackCountFilter.minTrackNPix = cms.int32(3)  # npix > 2 means >= 3
process.trackCountFilter.minNSelectedTracks = cms.int32(2)
process.trackCountFilter.maxNSelectedTracks = cms.int32(5)

process.eventFilter_HM = cms.Sequence(process.hltFilter * process.trackCountFilter)
process.eventFilter_HM_step = cms.Path(process.eventFilter_HM)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# ============================================================================
# ChiC Producers - All Modes Combined
# ============================================================================

# Common track quality cuts for all producers
COMMON_TRACK_CUTS = {
    'minTrackPt': cms.double(0.1),
    'maxTrackEta': cms.double(2.4),
    'maxTrackNormalizedChi2': cms.double(10.0),
    'minTrackNHits': cms.int32(5),  # nhits > 4
    'minTrackNPix': cms.int32(3),  # npix > 2
}

# --- 1. ChiC → 4π (ChiCFourTrackProducer) ---
from VertexCompositeAnalysis.VertexCompositeProducer.chiCTo4Pi_cfi import ChiCTo4Pi as _ChiCTo4Pi
process.ChiCTo4Pi = _ChiCTo4Pi.clone()
process.ChiCTo4Pi.trackCollection = cms.InputTag('generalTracks')
process.ChiCTo4Pi.daughterMasses = cms.vdouble(0.13957018, 0.13957018, 0.13957018, 0.13957018)
process.ChiCTo4Pi.daughterMassSigmas = cms.vdouble(0.0013957018, 0.0013957018, 0.0013957018, 0.0013957018)
for key, value in COMMON_TRACK_CUTS.items():
    setattr(process.ChiCTo4Pi, key, value)
process.ChiCTo4Pi.useVertexFitting = cms.bool(True)  # Enable vertex fitting
process.ChiCTo4Pi.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.ChiCTo4Pi.beamSpot = cms.InputTag('offlineBeamSpot')
process.ChiCTo4Pi.minCandidatePt = cms.double(0.0)
process.ChiCTo4Pi.minAcoplanarity = cms.double(0.0)
process.ChiCTo4Pi.maxSphericity = cms.double(1.0)
process.ChiCTo4Pi.maxCandidateAbsEta = cms.double(100.0)
process.ChiCTo4Pi.storeEventShape = cms.bool(True)
process.ChiCTo4Pi.applyMassWindow = cms.bool(False)  # No mass cut - store all candidates
process.ChiCTo4Pi.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC'),
        pdgId=cms.int32(445),
        mass=cms.double(3.35),  # Center of 2.5-4.2 range
        massWindow=cms.double(0.85)  # Covers 2.5-4.2 GeV (not used when applyMassWindow=False)
    )
)

# --- 2. Di-Kaon Producer (for ChiC → K+K- K+K-) ---
# DiKaonProducer reconstructs phi(1020) → K+K- resonances
# It requires: 2 tracks with opposite charge, kaon mass hypothesis, vertex fit
# Mass cuts: phi mass window (phiMass ± phiMassCut) AND mKK in [mKKCutMin, mKKCutMax]
from VertexCompositeAnalysis.VertexCompositeProducer.diKaonCandidates_cfi import DiKaonProducer as _DiKaonProducer
process.DiKaonProducer = _DiKaonProducer.clone()
process.DiKaonProducer.trackRecoAlgorithm = cms.InputTag('generalTracks')  # Will be replaced to unpackedTracksAndVertices by changeToMiniAOD
process.DiKaonProducer.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')  # Will be replaced to unpackedTracksAndVertices by changeToMiniAOD
process.DiKaonProducer.tkChi2Cut = cms.double(10.0)  # Track normalized chi2 < 10
process.DiKaonProducer.tkNhitsCut = cms.int32(5)  # Track nhits >= 5
process.DiKaonProducer.tkPtCut = cms.double(0.1)  # Track pt > 0.1 GeV
process.DiKaonProducer.tkEtaCut = cms.double(2.4)  # Track |eta| < 2.4
process.DiKaonProducer.phiMassCut = cms.double(1.0)  # Phi mass window: 1.019 ± 1.0 GeV (very loose, no tight cut)
process.DiKaonProducer.mKKCutMin = cms.double(0.2)  # K+K- invariant mass > 0.2 GeV (very loose)
process.DiKaonProducer.mKKCutMax = cms.double(2.0)  # K+K- invariant mass < 2.0 GeV (very loose)
# Loosen impact parameter cuts for debugging
process.DiKaonProducer.dauTransImpactSigCut = cms.double(0.0)  # Track transverse impact sig > 0 (very loose)
process.DiKaonProducer.dauLongImpactSigCut = cms.double(0.0)  # Track longitudinal impact sig > 0 (very loose)
process.DiKaonProducer.vtxSignificance3DCut = cms.double(0.0)  # 3D vertex significance > 0 (very loose)
process.DiKaonProducer.tkDCACut = cms.double(10.0)  # Track DCA < 10 cm (very loose)

# ChiC from Di-Kaons
from VertexCompositeAnalysis.VertexCompositeProducer.chiCFromDiKaons_cfi import ChiCFromDiKaons as _ChiCFromDiKaons
process.ChiCFromDiKaons = _ChiCFromDiKaons.clone()
process.ChiCFromDiKaons.resonanceCollection = cms.InputTag('DiKaonProducer', 'DiKaon')
process.ChiCFromDiKaons.applyMassWindow = cms.bool(False)  # No mass cut - store all candidates
process.ChiCFromDiKaons.requireUniqueTracks = cms.bool(True)
process.ChiCFromDiKaons.useVertexFitting = cms.bool(True)  # Enable vertex fitting
process.ChiCFromDiKaons.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.ChiCFromDiKaons.beamSpot = cms.InputTag('offlineBeamSpot')
process.ChiCFromDiKaons.resonanceMassSigmas = cms.vdouble(0.00493677, 0.00493677, 0.00493677, 0.00493677)  # Kaon mass uncertainty
process.ChiCFromDiKaons.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC'),
        pdgId=cms.int32(445),
        mass=cms.double(3.35),  # Center of 2.5-4.2 range
        massWindow=cms.double(0.85)  # Covers 2.5-4.2 GeV (not used when applyMassWindow=False)
    )
)

# --- 3. Kshort Producer (for ChiC → Ks Ks) ---
# KshortProducer reconstructs K0s → π+π- resonances
# It requires: 2 tracks with opposite charge, pion mass hypothesis, vertex fit
# Mass cuts: K0s mass window (0.497 ± kShortMassCut) AND π+π- mass in [mPiPiCutMin, mPiPiCutMax]
from VertexCompositeAnalysis.VertexCompositeProducer.kshortCandidates_cfi import KshortProducer as _KshortProducer
process.KshortProducer = _KshortProducer.clone()
process.KshortProducer.trackRecoAlgorithm = cms.InputTag('generalTracks')  # Will be replaced to unpackedTracksAndVertices by changeToMiniAOD
process.KshortProducer.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')  # Will be replaced to unpackedTracksAndVertices by changeToMiniAOD
process.KshortProducer.tkChi2Cut = cms.double(10.0)  # Track normalized chi2 < 10
process.KshortProducer.tkNhitsCut = cms.int32(5)  # Track nhits >= 5
process.KshortProducer.tkPtCut = cms.double(0.1)  # Track pt > 0.1 GeV
process.KshortProducer.tkEtaCut = cms.double(2.4)  # Track |eta| < 2.4
process.KshortProducer.tkDCACut = cms.double(1.0)  # Track DCA < 1.0 cm
process.KshortProducer.mPiPiCutMin = cms.double(0.0)  # π+π- invariant mass > 0.0 GeV
process.KshortProducer.mPiPiCutMax = cms.double(0.6)  # π+π- invariant mass < 0.6 GeV
process.KshortProducer.kShortMassCut = cms.double(0.030)  # K0s mass window: 0.497 ± 0.030 GeV (0.467-0.527 GeV)

# ChiC from Kshorts
from VertexCompositeAnalysis.VertexCompositeProducer.chiCFromKshorts_cfi import ChiCFromKshorts as _ChiCFromKshorts
process.ChiCFromKshorts = _ChiCFromKshorts.clone()
process.ChiCFromKshorts.resonanceCollection = cms.InputTag('KshortProducer', 'Kshort')
process.ChiCFromKshorts.applyMassWindow = cms.bool(False)  # No mass cut - store all candidates
process.ChiCFromKshorts.requireUniqueTracks = cms.bool(True)
process.ChiCFromKshorts.useVertexFitting = cms.bool(True)  # Enable vertex fitting
process.ChiCFromKshorts.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.ChiCFromKshorts.beamSpot = cms.InputTag('offlineBeamSpot')
process.ChiCFromKshorts.resonanceMassSigmas = cms.vdouble(0.0013957018, 0.0013957018, 0.0013957018, 0.0013957018)  # Pion mass uncertainty
process.ChiCFromKshorts.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC'),
        pdgId=cms.int32(445),
        mass=cms.double(3.35),  # Center of 2.5-4.2 range
        massWindow=cms.double(0.85)  # Covers 2.5-4.2 GeV (not used when applyMassWindow=False)
    )
)

# --- 4. Eta_c → p+ p- (ChiCTrackPairProducer) ---
from VertexCompositeAnalysis.VertexCompositeProducer.etaCToPP_cfi import EtaCToPP as _EtaCToPP
process.EtaCToPP = _EtaCToPP.clone()
process.EtaCToPP.trackCollection = cms.InputTag('generalTracks')
process.EtaCToPP.daughterMass = cms.double(0.938272013)  # Proton mass
process.EtaCToPP.daughterMassSigma = cms.double(0.00938272013)  # 1% uncertainty
for key, value in COMMON_TRACK_CUTS.items():
    setattr(process.EtaCToPP, key, value)
process.EtaCToPP.minPairPt = cms.double(0.0)
process.EtaCToPP.applyMassWindow = cms.bool(False)  # No mass cut - store all candidates
process.EtaCToPP.requiredChargeProduct = cms.int32(-1)  # Opposite charge (p+ p-)
process.EtaCToPP.useVertexFitting = cms.bool(True)  # Enable vertex fitting
process.EtaCToPP.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.EtaCToPP.beamSpot = cms.InputTag('offlineBeamSpot')
process.EtaCToPP.states = cms.VPSet(
    cms.PSet(
        name=cms.string('EtaC'),
        pdgId=cms.int32(441),
        mass=cms.double(3.35),  # Center of 2.5-4.2 range
        massWindow=cms.double(0.85)  # Covers 2.5-4.2 GeV (not used when applyMassWindow=False)
    )
)

# ============================================================================
# Flat Ntuplizer - All Modes
# ============================================================================
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCFlatNtuplizer_cfi import ChiCFlatNtuplizer as _ChiCFlatNtuplizer
process.ChiCFlatNtuplizer = _ChiCFlatNtuplizer.clone(
    treeName=cms.untracked.string('ChiCUltimateNtuple'),
    # MINIAOD uses offlineSlimmedPrimaryVertices - changeToMiniAOD will replace offlinePrimaryVertices 
    # with unpackedTracksAndVertices for producers, but ntuplizer should use the unpacked collection
    # After changeToMiniAOD, offlinePrimaryVertices points to unpackedTracksAndVertices which produces vertices
    primaryVertices=cms.InputTag('offlinePrimaryVertices'),  # Will be replaced by changeToMiniAOD to unpackedTracksAndVertices
    sources=cms.VPSet(
        cms.PSet(
            name=cms.string('ChiC_4Pi'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCTo4Pi', 'ChiC')
        ),
        cms.PSet(
            name=cms.string('ChiC_DiKaon'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCFromDiKaons', 'ChiC')
        ),
        cms.PSet(
            name=cms.string('ChiC_KsKs'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCFromKshorts', 'ChiC')
        ),
        cms.PSet(
            name=cms.string('EtaC_PP'),
            pdgId=cms.int32(441),
            collection=cms.InputTag('EtaCToPP', 'EtaC')
        )
    )
)

# ============================================================================
# Analysis Paths
# ============================================================================
process.chicUltimate_step = cms.Path(
    process.eventFilter_HM *
    process.ChiCTo4Pi *
    process.DiKaonProducer * process.ChiCFromDiKaons *
    process.KshortProducer * process.ChiCFromKshorts *
    process.EtaCToPP
)

process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff')
process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string('chic_ultimate_tree.root'),
)

process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
process.eventinfoana.triggerPathNames = cms.untracked.vstring('HLT_*')
process.eventinfoana.eventFilterNames = cms.untracked.vstring(
    'Flag_colEvtSel',
    'Flag_hfCoincFilter',
    'Flag_primaryVertexFilter',
)
process.eventinfoana.triggerFilterNames = cms.untracked.vstring()
process.eventinfoana.stageL1Trigger = cms.uint32(2)
process.pevt = cms.EndPath(process.eventinfoana)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.chicUltimate_step,
    process.pevt,
)

process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(
    process.eventFilter_HM * process.primaryVertexFilter * process.clusterCompatibilityFilter
)

eventFilterPaths = [process.Flag_colEvtSel, process.Flag_primaryVertexFilter]
for P in eventFilterPaths:
    process.schedule.insert(0, P)

changeToMiniAOD(process)
# Use unpackedTracksAndVertices for consistency with producers
# The unpacker (TrackAndVertexUnpacker) produces:
#   - reco::TrackCollection() [default instance]
#   - reco::VertexCollection() [default instance] - primary vertices with unpacked track associations
#   - reco::VertexCollection("secondary") - secondary vertices
# Using unpackedTracksAndVertices ensures vertices have proper track associations from unpacked tracks
process.ChiCFlatNtuplizer.primaryVertices = cms.InputTag('unpackedTracksAndVertices')
process.options.numberOfThreads = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# Enable Info level messages for debugging
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(1000))
process.MessageLogger.cerr.DiKaonProducer = cms.untracked.PSet(limit=cms.untracked.int32(1000))
process.MessageLogger.cerr.DiKaonFitter = cms.untracked.PSet(limit=cms.untracked.int32(1000))
process.MessageLogger.cerr.V0Fitter = cms.untracked.PSet(limit=cms.untracked.int32(1000))

# Output path for flat ntuple
process.chicNtuple = cms.EndPath(process.ChiCFlatNtuplizer)
process.schedule.append(process.chicNtuple)
