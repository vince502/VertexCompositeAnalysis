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
        description='ChiC Analysis Configuration for Condor',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cmsRun pp2024_ChiC_MB_cfg_v1_ForCondor.py inputFiles=file1.root,file2.root outputFile=output.root
  cmsRun pp2024_ChiC_MB_cfg_v1_ForCondor.py inputFileList=filelist.txt outputFile=output.root storageSite=cern
  cmsRun pp2024_ChiC_MB_cfg_v1_ForCondor.py inputFileList=filelist.txt outputFile=output.root storageSite=fnal
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
    #'infn': 'root://xrootd-se37-vanderbilt.sites.opensciencegrid.org/',
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
'file:04e18742-3308-45a5-b0d6-560741bec33f.root',
#        '/store/hidata/OORun2025/IonPhysics0/MINIAOD/PromptReco-v1/000/394/075/00000/09db905b-c8ac-4e9e-9d6d-2be7f844a12b.root'
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
print("ChiC Analysis Configuration")
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
    # Optional: specify storage account/site
    # This helps CMSSW route file access to the correct storage element
    # StorageAccount = cms.untracked.string('site://T2_CH_CERN'),  # Example for CERN
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('150X_dataRun3_Prompt_v1')

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_*',
]

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')
process.colEvtSel = cms.Sequence()

process.eventFilter_HM = cms.Sequence(process.hltFilter)
process.eventFilter_HM_step = cms.Path(process.eventFilter_HM)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# --- Resonance producers ---
from VertexCompositeAnalysis.VertexCompositeProducer.chiCTo4Pi_cfi import ChiCTo4Pi as _ChiCTo4Pi
from VertexCompositeAnalysis.VertexCompositeProducer.chiC2To4K_cfi import ChiC2To4K as _ChiC2To4K
from VertexCompositeAnalysis.VertexCompositeProducer.chiCTo2Ka_cfi import ChiCTo2Ka as _ChiCTo2Ka
from VertexCompositeAnalysis.VertexCompositeProducer.diKaonCandidates_cfi import DiKaonProducer as _DiKaonProducer
from VertexCompositeAnalysis.VertexCompositeProducer.chiCFromDiKaons_cfi import ChiCFromDiKaons as _ChiCFromDiKaons
from VertexCompositeAnalysis.VertexCompositeProducer.kshortCandidates_cfi import KshortProducer as _KshortProducer
from VertexCompositeAnalysis.VertexCompositeProducer.chiCFromKshorts_cfi import ChiCFromKshorts as _ChiCFromKshorts
from VertexCompositeAnalysis.VertexCompositeAnalyzer.chiCNtuplizer_cfi import ChiCNtuplizer as _ChiCNtuplizer

process.ChiCTo4Pi = _ChiCTo4Pi.clone()
# Expose key ChiCTo4Pi selections
process.ChiCTo4Pi.minTrackPt = cms.double(2)
process.ChiCTo4Pi.maxTrackEta = cms.double(1.6)
process.ChiCTo4Pi.maxTrackNormalizedChi2 = cms.double(10.0)
process.ChiCTo4Pi.minTrackNHits = cms.int32(0)
process.ChiCTo4Pi.minCandidatePt = cms.double(5)
process.ChiCTo4Pi.minAcoplanarity = cms.double(0.6)
process.ChiCTo4Pi.maxSphericity = cms.double(0.35)
process.ChiCTo4Pi.maxCandidateAbsEta = cms.double(2.0)
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
process.ChiC2To4K = _ChiC2To4K.clone()
# Expose key ChiC2To4K selections for easy tweaking in this config
process.ChiC2To4K.minTrackPt = cms.double(2)
process.ChiC2To4K.maxTrackEta = cms.double(2.4)
process.ChiC2To4K.maxTrackNormalizedChi2 = cms.double(10.0)
process.ChiC2To4K.minTrackNHits = cms.int32(0)
process.ChiC2To4K.minCandidatePt = cms.double(0.0)
process.ChiC2To4K.minAcoplanarity = cms.double(0.0)
process.ChiC2To4K.maxSphericity = cms.double(1.5)
process.ChiC2To4K.maxCandidateAbsEta = cms.double(2.4)
process.ChiC2To4K.storeEventShape = cms.bool(False)
process.ChiC2To4K.applyMassWindow = cms.bool(True)
process.ChiC2To4K.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC0'),
        pdgId=cms.int32(10441),
        mass=cms.double(3.4147),
        massWindow=cms.double(0.190)
    ),
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.190)
    )
)
process.ChiCTo2Ka = _ChiCTo2Ka.clone()
# Expose key ChiCTo2Ka selections
process.ChiCTo2Ka.minTrackPt = cms.double(2.0)
process.ChiCTo2Ka.maxTrackEta = cms.double(2.4)
process.ChiCTo2Ka.maxTrackNormalizedChi2 = cms.double(10.0)
process.ChiCTo2Ka.minTrackNHits = cms.int32(0)
process.ChiCTo2Ka.minPairPt = cms.double(5)
process.ChiCTo2Ka.applyMassWindow = cms.bool(True)
process.ChiCTo2Ka.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC0'),
        pdgId=cms.int32(10441),
        mass=cms.double(3.4147),
        massWindow=cms.double(0.180)
    ),
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.180)
    )
)

process.KshortProducer = _KshortProducer.clone()
# Expose key KshortProducer selections
process.KshortProducer.trackRecoAlgorithm = cms.InputTag('generalTracks')
process.KshortProducer.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.KshortProducer.tkChi2Cut = cms.double(10.0)
process.KshortProducer.tkNhitsCut = cms.int32(0)
process.KshortProducer.tkPtCut = cms.double(1.0)
process.KshortProducer.tkDCACut = cms.double(1.0)
process.KshortProducer.mPiPiCutMin = cms.double(0.0)
process.KshortProducer.mPiPiCutMax = cms.double(0.6)
process.KshortProducer.kShortMassCut = cms.double(0.030)
process.ChiC2FromKshorts = _ChiCFromKshorts.clone(
    states=cms.VPSet(
        cms.PSet(
            name=cms.string('ChiC2'),
            pdgId=cms.int32(445),
            mass=cms.double(3.5562),
            massWindow=cms.double(0.150)
        )
    )
)
process.ChiC2FromKshortsSequence = cms.Sequence(process.KshortProducer * process.ChiC2FromKshorts)

# Di-kaon pair combination (wide mass range)
process.DiKaonWideMass = _DiKaonProducer.clone(
    phiMassCut=cms.double(0.7),
    mKKCutMin=cms.double(0.4),
    mKKCutMax=cms.double(1.5)
)
# Expose di-kaon mass window knobs (0.4-1.5 GeV default)
process.DiKaonWideMass.phiMassCut = cms.double(0.7)
process.DiKaonWideMass.mKKCutMin = cms.double(0.4)
process.DiKaonWideMass.mKKCutMax = cms.double(1.5)

process.ChiC2FromDiKaonPairs = _ChiCFromDiKaons.clone(
    resonanceCollection=cms.InputTag('DiKaonWideMass', 'DiKaon')
)
process.ChiC2FromDiKaonPairs.applyMassWindow = cms.bool(True)
process.ChiC2FromDiKaonPairs.requireUniqueTracks = cms.bool(True)
process.ChiC2FromDiKaonPairs.states = cms.VPSet(
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.160)
    )
)
process.ChiC2FromDiKaonPairsSequence = cms.Sequence(process.DiKaonWideMass * process.ChiC2FromDiKaonPairs)

# Unified ChiC ntuple writer
# cand_type mapping (0-based index in sources list):
#   cand_type = 0: ChiC0 → 4π (ChiC0_4Pi)
#   cand_type = 1: ChiC2 → 4π (ChiC2_4Pi)
#   cand_type = 2: ChiC0 → 4K (ChiC0_4K)
#   cand_type = 3: ChiC2 → 4K (ChiC2_4K)
#   cand_type = 4: ChiC0 → 2Ka (ChiC0_2Ka)
#   cand_type = 5: ChiC2 → 2Ka (ChiC2_2Ka)
#   cand_type = 6: ChiC2 → Kshort Kshort (ChiC2_KshortKshort)
#   cand_type = 7: ChiC2 → KK (ChiC2_KK)
process.ChiCNtuplizer = _ChiCNtuplizer.clone(
    treeName=cms.untracked.string('ChiCNtuple'),
    storeDaughterInfo=cms.untracked.bool(True),
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
        ),
        cms.PSet(  # cand_type = 2: ChiC0 → 4K
            name=cms.string('ChiC0_4K'),
            pdgId=cms.int32(10441),
            collection=cms.InputTag('ChiC2To4K', 'ChiC0')
        ),
        cms.PSet(  # cand_type = 3: ChiC2 → 4K
            name=cms.string('ChiC2_4K'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiC2To4K', 'ChiC2')
        ),
        cms.PSet(  # cand_type = 4: ChiC0 → 2Ka
            name=cms.string('ChiC0_2Ka'),
            pdgId=cms.int32(10441),
            collection=cms.InputTag('ChiCTo2Ka', 'ChiC0')
        ),
        cms.PSet(  # cand_type = 5: ChiC2 → 2Ka
            name=cms.string('ChiC2_2Ka'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCTo2Ka', 'ChiC2')
        ),
        cms.PSet(  # cand_type = 6: ChiC2 → Kshort Kshort
            name=cms.string('ChiC2_KshortKshort'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiC2FromKshorts', 'ChiC2')
        ),
        cms.PSet(  # cand_type = 7: ChiC2 → KK
            name=cms.string('ChiC2_KK'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiC2FromDiKaonPairs', 'ChiC2')
        )
    )
)

# --- Analysis paths ---
process.chic4Pi_step = cms.Path(
    process.eventFilter_HM * process.ChiCTo4Pi
)

process.chic4K_step = cms.Path(
    process.eventFilter_HM * process.ChiC2To4K
)

process.chic2K_step = cms.Path(
    process.eventFilter_HM * process.ChiC2FromKshortsSequence
)

process.chic2KK_step = cms.Path(
    process.eventFilter_HM * process.ChiC2FromDiKaonPairsSequence
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
    process.chic4Pi_step,
    process.chic4K_step,
    process.chic2K_step,
    process.chic2KK_step,
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
process.chicNtuple = cms.EndPath(process.ChiCNtuplizer)
#process.outpathcustom = cms.EndPath(process.outCustom)
#process.schedule.append(process.outpathcustom)
process.schedule.append(process.chicNtuple)
