"""
Configuration for OniaPhotonProducer
Producer for Onia + Photon (conversion) candidates

Reconstructs:
  - chi_c -> J/psi + gamma  (deltaM ~ 0.3-0.5 GeV)
  - chi_b -> Upsilon + gamma (deltaM ~ 0.4-0.8 GeV)
  - D*0 -> D0 + gamma (deltaM ~ 0.14 GeV, rare radiative decay)

Input:
  - oniaSrc: J/psi, Upsilon, or D0 candidates (pat::CompositeCandidate)
  - conversionSrc: Conversion photons from OniaPhotonConversionProducer

Output:
  - pat::CompositeCandidateCollection with chi candidates
  - Optionally: kinematically refitted candidates (with mass constraint)
"""

import FWCore.ParameterSet.Config as cms

# Base configuration for chi_c -> J/psi + gamma
chiCandidates = cms.EDProducer(
    "OniaPhotonProducer",
    # Input collections
    oniaSrc=cms.InputTag("onia2MuMuPAT"),
    conversionSrc=cms.InputTag("oniaPhotonCandidates", "conversions"),
    pvSrc=cms.InputTag("offlinePrimaryVertices"),
    beamSpotSrc=cms.InputTag("offlineBeamSpot"),
    # Conversion cuts
    minConvPt=cms.double(0.0),
    minConvRho=cms.double(0.0),
    rejectPi0=cms.bool(False),
    # Onia cuts
    minOniaPt=cms.double(0.0),
    maxOniaAbsY=cms.double(2.4),
    requireTriggerMatch=cms.bool(False),
    # Matching cuts
    dzMax=cms.double(1.0),
    deltaMassRange=cms.vdouble(0.0, 2.0),
    # Combined candidate cuts
    minCandPt=cms.double(0.0),
    maxCandAbsY=cms.double(2.4),
    # Kinematic refit options
    doKinematicRefit=cms.bool(False),
    constraintMass=cms.double(3.0969),
    # Parent type: "jpsi", "upsilon1S", "upsilon2S", "upsilon3S", "d0"
    parentType=cms.string("jpsi"),
)

# chi_c with kinematic refit
chiCandidatesRefit = chiCandidates.clone(
    doKinematicRefit=cms.bool(True),
)

# chi_b -> Upsilon(1S) + gamma
chibCandidates1S = chiCandidates.clone(
    constraintMass=cms.double(9.46030),
    deltaMassRange=cms.vdouble(0.0, 2.0),
    parentType=cms.string("upsilon1S"),
)

# chi_b -> Upsilon(2S) + gamma
chibCandidates2S = chiCandidates.clone(
    constraintMass=cms.double(10.02326),
    deltaMassRange=cms.vdouble(0.0, 2.0),
    parentType=cms.string("upsilon2S"),
)

# chi_b -> Upsilon(3S) + gamma
chibCandidates3S = chiCandidates.clone(
    constraintMass=cms.double(10.35520),
    deltaMassRange=cms.vdouble(0.0, 2.0),
    parentType=cms.string("upsilon3S"),
)

# D*0 -> D0 + gamma (radiative decay)
dstar0Candidates = chiCandidates.clone(
    oniaSrc=cms.InputTag("generalD0Candidates"),
    constraintMass=cms.double(1.86484),
    deltaMassRange=cms.vdouble(0.0, 0.5),
    parentType=cms.string("d0"),
    doKinematicRefit=cms.bool(False),
)

# Heavy-ion specific configuration (tighter cuts)
chiCandidatesHI = chiCandidates.clone(
    minOniaPt=cms.double(6.5),
    minCandPt=cms.double(8.0),
    dzMax=cms.double(0.5),
    deltaMassRange=cms.vdouble(0.0, 1.0),
)

# Loose configuration for studies
chiCandidatesLoose = chiCandidates.clone(
    dzMax=cms.double(2.0),
    deltaMassRange=cms.vdouble(0.0, 3.0),
)


def configureChiCProduction(process, isMC=False, isHI=False, doRefit=True):
    """
    Configure chi_c production in a CMSSW process.
    """
    _ = isMC
    baseCfg = chiCandidatesHI if isHI else chiCandidates
    process.chiProducer = baseCfg.clone(
        doKinematicRefit=cms.bool(doRefit),
    )
    return cms.Sequence(process.chiProducer)
