import FWCore.ParameterSet.Config as cms

chiCJpsis = cms.EDProducer(
    "ChiCJpsiProducer",
    muons=cms.InputTag("slimmedMuons"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot=cms.InputTag("offlineBeamSpot"),
    minMuonPt=cms.double(2.0),
    maxMuonAbsEta=cms.double(2.4),
    minMass=cms.double(2.7),
    maxMass=cms.double(3.5),
    minVertexProbability=cms.double(0.005),
)

# Converts legacy Onia2MuMu candidates into explicit vertex-fitted states with
# the full (x,y,z,px,py,pz,m) covariance. It is also the sole source of the
# prompt common-origin hypothesis used by the one-prong branch.
chiCOniaStates = cms.EDProducer(
    "ChiCOniaStateProducer",
    jpsis=cms.InputTag("chiCJpsis"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    useEmbeddedPVIndex=cms.bool(True),
    useOutputIndices=cms.bool(False),
    oniaPDGs=cms.vint32(),
    massRanges=cms.vdouble(),
    nominalMasses=cms.vdouble(),
    oppositeSign=cms.bool(False),
)

chiCConversionLegs = cms.EDProducer(
    "ChiCConversionLegProducer",
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    minPt=cms.double(0.2),
    maxAbsEta=cms.double(2.5),
    requireHighPurity=cms.bool(False),
)

chiCSingleLegGeometries = cms.EDProducer(
    "ChiCSingleLegGeometryProducer",
    jpsis=cms.InputTag("chiCOniaStates"),
    legs=cms.InputTag("chiCConversionLegs"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    maxDeltaRToJpsi=cms.double(1.5),
    minConversionRadius=cms.double(1.0),
    maxConversionRadius=cms.double(110.0),
    maxAbsConversionZ=cms.double(280.0),
    maxBackwardPath=cms.double(100.0),
    nGeometrySteps=cms.uint32(64),
    minDirectionCosine=cms.double(0.90),
    maxRoots=cms.uint32(3),
)

chiCOneProngCandidates = cms.EDProducer(
    "ChiCOneProngCandidateProducer",
    geometries=cms.InputTag("chiCSingleLegGeometries"),
    jpsis=cms.InputTag("chiCOniaStates"),
    legs=cms.InputTag("chiCConversionLegs"),
)

chiCConventionalPhotonStates = cms.EDProducer(
    "ChiCConventionalPhotonStateProducer",
    candidates=cms.InputTag("chiCConventionalCandidates"),
)
