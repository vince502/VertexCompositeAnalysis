import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dStar5panalyzer_tree_cfi import *

# Wrong-sign configuration expects matching producer/module to be defined upstream.
dStar5pana_wrongsign = dStar5pana.clone(
    CompositeCollection=cms.untracked.InputTag('generalDStar5PCandidatesWrongSign:DStar5P'),
    MVACollection=cms.InputTag('generalDStar5PCandidatesWrongSign:MVAValuesDStar5P')
)

dStar5pana_mc_wrongsign = dStar5pana_mc.clone(
    CompositeCollection=cms.untracked.InputTag('generalDStar5PCandidatesWrongSign:DStar5P'),
    MVACollection=cms.InputTag('generalDStar5PCandidatesWrongSign:MVAValuesDStar5P')
)

# TOF matching variants

dStar5pana_tof = dStar5pana.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)

dStar5pana_tof_wrongsign = dStar5pana_wrongsign.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)

dStar5pana_tof_mc = dStar5pana_mc.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)

dStar5pana_tof_mc_wrongsign = dStar5pana_mc_wrongsign.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)
