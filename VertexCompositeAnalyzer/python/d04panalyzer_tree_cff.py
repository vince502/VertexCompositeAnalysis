import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.d04panalyzer_tree_cfi import *

# Wrong-sign combination assumes dedicated producer/job provides the collection
d04pana_wrongsign = d04pana.clone(
    CompositeCollection=cms.untracked.InputTag('generalD04PCandidatesWrongSign:D04P'),
    MVACollection=cms.InputTag('generalD04PCandidatesWrongSign:MVAValuesD04P')
)

d04pana_mc_wrongsign = d04pana_mc.clone(
    CompositeCollection=cms.untracked.InputTag('generalD04PCandidatesWrongSign:D04P'),
    MVACollection=cms.InputTag('generalD04PCandidatesWrongSign:MVAValuesD04P')
)

# Variants enabling TOF matching

d04pana_tof = d04pana.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)

d04pana_tof_wrongsign = d04pana_wrongsign.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)

d04pana_tof_mc = d04pana_mc.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)

d04pana_tof_mc_wrongsign = d04pana_mc_wrongsign.clone(
    doGenMatchingTOF=cms.untracked.bool(True)
)
