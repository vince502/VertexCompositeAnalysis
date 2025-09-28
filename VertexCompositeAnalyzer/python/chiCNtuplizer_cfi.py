import FWCore.ParameterSet.Config as cms

ChiCNtuplizer = cms.EDAnalyzer(
    'ChiCNtuplizer',
    treeName=cms.untracked.string('ChiCNtuple'),
    storeDaughterInfo=cms.untracked.bool(True),
    sources=cms.VPSet(
        # To be configured in the job file.
    )
)
