import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import sys
import os
import argparse

# Parse command line arguments using argparse
# CMSSW passes arguments in format: key=value key2=value2
# We need to parse these before argparse can handle them
def parse_cmssw_args():
    """Parse CMSSW-style arguments (key=value) into argparse format"""
    parser = argparse.ArgumentParser(
        description='UPC2025 ChiC Analysis Configuration for Condor',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cmsRun UPC2025_ChiC_MB_cfg_v1_ForCondor.py inputFiles=file1.root,file2.root outputFile=output.root
  cmsRun UPC2025_ChiC_MB_cfg_v1_ForCondor.py inputFileList=filelist.txt outputFile=output.root storageSite=cern
  cmsRun UPC2025_ChiC_MB_cfg_v1_ForCondor.py inputFileList=filelist.txt outputFile=output.root storageSite=fnal
        """
    )
    
    parser.add_argument('--inputFiles', type=str, default=None,
                       help='Comma-separated list of input ROOT files')
    parser.add_argument('--inputFileList', type=str, default=None,
                       help='Path to text file containing list of input files (one per line)')
    parser.add_argument('--outputFile', type=str, default='chic_combination_tree.root',
                       help='Output ROOT file path (default: chic_combination_tree.root)')
    parser.add_argument('--storageSite', type=str, default='cern',
                       help='Storage site/redirector: cern (default), fnal, local, or custom xrootd redirector (e.g., root://hostname/)')
    
    # Parse CMSSW-style arguments (key=value format)
    # Convert to standard --key=value format for argparse
    cmssw_args = []
    for arg in sys.argv[1:]:
        if '=' in arg and not arg.startswith('--'):
            # Convert key=value to --key=value
            key, value = arg.split('=', 1)
            cmssw_args.append(f'--{key}={value}')
        else:
            cmssw_args.append(arg)
    
    args = parser.parse_args(cmssw_args)
    return args

# Parse arguments
args = parse_cmssw_args()

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

# Get storage site/redirector from parsed arguments
storageSite = args.storageSite.lower()

# Define xrootd redirectors for different sites
xrootdRedirectors = {
    'cern': 'root://eoscms.cern.ch/',
    'fnal': 'root://cmsxrootd.fnal.gov/',
    'infn': 'root://xrootd-cms.infn.it/',
    'local': '',  # No prefix for local files
}

# Get redirector prefix
if storageSite in xrootdRedirectors:
    redirectorPrefix = xrootdRedirectors[storageSite]
else:
    # Custom redirector provided (not in predefined list)
    redirectorPrefix = storageSite if storageSite.startswith('root://') else 'root://' + storageSite + '/'

# Get input files from parsed arguments
inputFiles = []

if args.inputFiles:
    # Comma-separated list provided
    inputFiles = [f.strip() for f in args.inputFiles.split(',') if f.strip()]
elif args.inputFileList:
    # File list provided
    fileListPath = args.inputFileList
    if os.path.exists(fileListPath):
        with open(fileListPath, 'r') as f:
            inputFiles = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    else:
        raise FileNotFoundError(f"Input file list not found: {fileListPath}")

# Default if no arguments provided
if not inputFiles:
    inputFiles = [
        '/store/hidata/HIRun2025A/HIForward0/MINIAOD/PromptReco-v1/000/399/540/00000/491ce449-be44-4fe8-a337-f85c90499ea9.root',
    ]

# Add redirector prefix to file paths if needed
# Only add prefix if file path starts with /store/ and redirector is not 'local'
processedFiles = []
for f in inputFiles:
    if f.startswith('/store/') and redirectorPrefix and storageSite != 'local':
        # Add xrootd redirector prefix
        processedFiles.append(redirectorPrefix + f)
    elif f.startswith('root://'):
        # Already has redirector, use as is
        processedFiles.append(f)
    elif f.startswith('file:'):
        # Local file, use as is
        processedFiles.append(f)
    else:
        # Assume it's a local path or already has redirector
        processedFiles.append(f)

inputFiles = processedFiles

# Get output file from parsed arguments
outputFile = args.outputFile

# Print configuration summary
print("=" * 80)
print("UPC2025 ChiC Analysis Configuration")
print("=" * 80)
print("Storage site: {0}".format(storageSite))
if redirectorPrefix:
    print("Xrootd redirector: {0}".format(redirectorPrefix))
print("Input files ({0}):".format(len(inputFiles)))
for i, f in enumerate(inputFiles[:5], 1):  # Print first 5
    print("  [{0}] {1}".format(i, f))
if len(inputFiles) > 5:
    print("  ... and {0} more files".format(len(inputFiles) - 5))
print("Output file: {0}".format(outputFile))
print("=" * 80)

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(inputFiles),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('151X_dataRun3_Prompt_v1')

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v16*',
]

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')
process.colEvtSel = cms.Sequence()

# Track count filter - configurable cut on number of selected tracks
from VertexCompositeAnalysis.VertexCompositeProducer.trackCountFilter_cfi import trackCountFilter
process.trackCountFilter = trackCountFilter.clone()
# Configure track selection cuts (matching ChiC producers)
process.trackCountFilter.trackCollection = cms.InputTag('generalTracks')
process.trackCountFilter.minTrackPt = cms.double(0.1)  # Match ChiCTo4Pi/ChiCTo2Ka
process.trackCountFilter.maxTrackEta = cms.double(2.4)  # Match ChiC2To4K/ChiCTo2Ka
process.trackCountFilter.maxTrackNormalizedChi2 = cms.double(10.0)
process.trackCountFilter.minTrackNHits = cms.int32(0)  # Match ChiCTo4Pi/ChiC2To4K/ChiCTo2Ka
# Configure track count cuts
process.trackCountFilter.minNSelectedTracks = cms.int32(4)  # Minimum tracks required (>=)
process.trackCountFilter.maxNSelectedTracks = cms.int32(6)  # Maximum tracks (<=)

process.eventFilter_HM = cms.Sequence(process.hltFilter * process.trackCountFilter)
process.eventFilter_HM_step = cms.Path(process.eventFilter_HM)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# --- D0 and B producers ---
process.load('VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff')
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(0)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.3)
process.generalD0CandidatesNew.tkEtaCut = cms.double(2.4)
process.generalD0CandidatesNew.tkChi2Cut = cms.double(10.0)
process.generalD0CandidatesNew.d0MassCut = cms.double(0.15)

process.load('VertexCompositeAnalysis.VertexCompositeProducer.generalBCandidates_cff')
process.generalBCandidatesNew = process.generalBCandidates.clone()
process.generalBCandidatesNew.d0RecoAlgorithm = cms.InputTag('generalD0CandidatesNew', 'D0')
process.generalBCandidatesNew.batTkPtCut = cms.double(0.3)
process.generalBCandidatesNew.batTkEtaCut = cms.double(2.4)
process.generalBCandidatesNew.batTkChi2Cut = cms.double(10.0)
process.generalBCandidatesNew.batTkNhitsCut = cms.int32(0)

# --- ChiC producers (4pi only) ---
from VertexCompositeAnalysis.VertexCompositeProducer.chiCTo4Pi_cfi import ChiCTo4Pi as _ChiCTo4Pi
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCFlatNtuplizer_cfi import ChiCFlatNtuplizer as _ChiCFlatNtuplizer

process.ChiCTo4Pi = _ChiCTo4Pi.clone()
# Expose key ChiCTo4Pi selections
process.ChiCTo4Pi.minTrackPt = cms.double(0.1)
process.ChiCTo4Pi.maxTrackEta = cms.double(2.4)
process.ChiCTo4Pi.maxTrackNormalizedChi2 = cms.double(999.0)
process.ChiCTo4Pi.minTrackNHits = cms.int32(0)
process.ChiCTo4Pi.minCandidatePt = cms.double(0)
process.ChiCTo4Pi.minAcoplanarity = cms.double(0.0)
process.ChiCTo4Pi.maxSphericity = cms.double(1.0)
process.ChiCTo4Pi.maxCandidateAbsEta = cms.double(2.4)
process.ChiCTo4Pi.storeEventShape = cms.bool(True)
process.ChiCTo4Pi.applyMassWindow = cms.bool(True)
process.ChiCTo4Pi.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC0'),
        pdgId=cms.int32(10441),
        mass=cms.double(3.4147),
        massWindow=cms.double(0.125)
    ),
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.125)
    )
)

# Unified ChiC flat ntuple writer (4pi only) - comprehensive information for offline analysis
# cand_type mapping (0-based index in sources list):
#   cand_type = 0: ChiC0 → 4π (ChiC0_4Pi)
#   cand_type = 1: ChiC2 → 4π (ChiC2_4Pi)
process.ChiCFlatNtuplizer = _ChiCFlatNtuplizer.clone(
    treeName=cms.untracked.string('ChiCFlatNtuple'),
    primaryVertices=cms.InputTag('offlinePrimaryVertices'),
    sources=cms.VPSet(
        cms.PSet(  # cand_type = 0: ChiC0 → 4π
            name=cms.string('ChiC0_4Pi'),
            pdgId=cms.int32(10441),
            collection=cms.InputTag('ChiCTo4Pi', 'ChiC0')
        ),
        cms.PSet(  # cand_type = 1: ChiC2 → 4π
            name=cms.string('ChiC2_4Pi'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCTo4Pi', 'ChiC2')
        )
    )
)

# --- Analysis paths ---
# D0 and B reconstruction
process.d0_step = cms.Path(
    process.eventFilter_HM * process.generalD0CandidatesNew
)

process.b_step = cms.Path(
    process.eventFilter_HM * process.generalD0CandidatesNew * process.generalBCandidatesNew
)

# ChiC to 4pi only
process.chic4Pi_step = cms.Path(
    process.eventFilter_HM * process.ChiCTo4Pi
)

process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff')
process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string(outputFile),
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
    #process.d0_step,
    #process.b_step,
    process.chic4Pi_step,
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
process.options.numberOfThreads = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.outCustom = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("myOutput.root"),
    outputCommands = cms.untracked.vstring("keep *_*_*_ANASKIM")  # Keep everything
)

# Output path
process.chicNtuple = cms.EndPath(process.ChiCFlatNtuplizer)
#process.outpathcustom = cms.EndPath(process.outCustom)
#process.schedule.append(process.outpathcustom)
process.schedule.append(process.chicNtuple)
