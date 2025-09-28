import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalV0Candidates_cfi import generalV0Candidates as _generalV0Candidates

DiKaonCandidates = _generalV0Candidates.clone(
    selectKshorts = False,
    selectLambdas = False,
    selectPhis = True,
    selectD0s = False,
    selectDSToKsKs = False,
    selectDSToPhiPis = False,
    selectDPMs = False,
    selectLambdaCToLamPis = False,
    selectLambdaCToKsPs = False,
    selectXis = False,
    selectOmegas = False,
    phiMassCut = 0.02
)

_DiKaonParams = {name: getattr(DiKaonCandidates, name) for name in DiKaonCandidates.parameterNames_()}

DiKaonProducer = cms.EDProducer(
    "DiKaonProducer",
    **_DiKaonParams
)
