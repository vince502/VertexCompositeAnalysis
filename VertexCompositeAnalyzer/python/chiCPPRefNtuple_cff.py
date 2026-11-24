import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCCombinedMiniAOD_cff import *
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCEventInfoNtuplizer_cfi import (
    chiCEventInfoNtuplizer,
)

# Shared 2024 pp-reference reconstruction selections. Both data and MC use
# the same reconstructed-object selection so their ntuples remain comparable.
chiCMinPtSlimmedMuons.cut = cms.string("pt > 2.0")
chiCMinPtSlimmedMuonCounter.minNumber = cms.uint32(2)

chiCSharedOniaSelectedMuons.cut = cms.string(
    'muonID("AllTrackerMuons")'
    ' && muonID("TMOneStationTight")'
    " && abs(innerTrack.dxy) < 0.3"
    " && abs(innerTrack.dz) < 20."
    " && innerTrack.hitPattern.trackerLayersWithMeasurement > 5"
    " && innerTrack.hitPattern.pixelLayersWithMeasurement > 0"
    " && pt > 2.0"
    " && abs(eta) <= 2.4"
    " && ("
    "      (abs(eta) < 1.0 && pt > 3.3)"
    "   || (abs(eta) >= 1.0 && abs(eta) < 1.3"
    "       && pt > 3.3 - 4.0 * (abs(eta) - 1.0))"
    "   || (abs(eta) >= 1.3 && abs(eta) < 1.7"
    "       && pt > 1.53 - 1.325 * (abs(eta) - 1.3))"
    "   || (abs(eta) >= 1.7 && abs(eta) <= 2.4 && pt > 1.0)"
    " )"
)

chiCSharedOniaJpsis.dimuonSelection = cms.string(
    "2.7 < mass && mass < 3.5 && pt < 50. && abs(y) < 2.4"
    " && charge == 0 && userFloat('vProb') > 0.01"
)

# Data uses the established fully reconstructed two-prong conversion branch.
chiCPPRefDataCovarianceNtuplizer = chiCCombinedCovarianceNtuplizer.clone(
    doOneProng=cms.bool(False),
    isMC=cms.bool(False),
)

chiCPPRefDataEventInfoNtuplizer = chiCEventInfoNtuplizer.clone(
    isMC=cms.bool(False),
)

chiCPPRefDataSequence = cms.Sequence(
    chiCMinPtSlimmedMuonFilterSequence
    * chiCTwoProngMiniAODJpsiFilteredSequence
)

# MC additionally runs the one-prong hypothesis and generator diagnostics.
chiCPPRefMCCovarianceNtuplizer = chiCCombinedCovarianceNtuplizer.clone(
    doOneProng=cms.bool(True),
    isMC=cms.bool(True),
)

chiCPPRefMCGenTruthNtuplizer = chiCGenTruthNtuplizer.clone()

chiCPPRefMCEventInfoNtuplizer = chiCEventInfoNtuplizer.clone(
    isMC=cms.bool(True),
)

chiCPPRefMCSequence = cms.Sequence(
    chiCMinPtSlimmedMuonFilterSequence
    * chiCSharedOniaMiniAODJpsiFilteredSequence
    * chiCConventionalMiniAODSequence
    * chiCOneProngMiniAODSequence
)
