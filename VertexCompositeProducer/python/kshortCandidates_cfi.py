import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalV0Candidates_cfi import generalV0Candidates as _generalV0Candidates

KshortCandidates = _generalV0Candidates.clone(
    selectKshorts = True,
    selectLambdas = False,
    selectPhis = False,
    selectD0s = False,
    selectDSToKsKs = False,
    selectDSToPhiPis = False,
    selectDPMs = False,
    selectLambdaCToLamPis = False,
    selectLambdaCToKsPs = False,
    selectXis = False,
    selectOmegas = False,
    kShortMassCut = 0.030
)

_KshortParams = {name: getattr(KshortCandidates, name) for name in KshortCandidates.parameterNames_()}

KshortProducer = cms.EDProducer(
    "KshortProducer",
    **_KshortParams
)
