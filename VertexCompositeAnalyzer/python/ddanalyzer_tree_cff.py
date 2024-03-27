import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.ddanalyzer_tree_cfi import *

ddana_wrongsign = ddana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalDStarCandidatesNewWrongSign:D0")
  #MVACollection = cms.InputTag("generalDStarCandidatesNewWrongSign:MVAValues")
)

#ddana_wrongsign_mc = ddana_mc.clone(
#  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
#  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
#)

ddana_tof = ddana.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

#ddana_tof_wrongsign = ddana_wrongsign.clone(
#  doGenMatchingTOF = cms.untracked.bool(True)
#)

ddana_tof_mc = ddana_mc.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

#ddana_tof_wrongsign_mc = ddana_wrongsign_mc.clone(
#  doGenMatchingTOF = cms.untracked.bool(True)
#)

