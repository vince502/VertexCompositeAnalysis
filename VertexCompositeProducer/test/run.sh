#!/bin/bash

# Script to run CMSSW job with a subset of files from file.txt
# Usage: ./run.sh <job_index> <files_per_job> <output_dir>

JOB_INDEX=$1
FILES_PER_JOB=${2:-20}
OUTPUT_DIR=${3:-/afs/cern.ch/work/s/soohwan/private/Analysis/DmesonAna/2025OxygenAnalysis/output}

# Calculate line range for this job
START_LINE=$((JOB_INDEX * FILES_PER_JOB + 1))
END_LINE=$((START_LINE + FILES_PER_JOB - 1))

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/pp2024_ChiC_MB_cfg_v1.py"
FILE_LIST="${SCRIPT_DIR}/file.txt"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Create temporary config file with selected files
TEMP_CONFIG="${SCRIPT_DIR}/pp2024_ChiC_MB_cfg_v${JOB_INDEX}.py"
OUTPUT_FILE="${OUTPUT_DIR}/chic_combination_tree_job${JOB_INDEX}.root"

# Extract files for this job
JOB_FILES=$(sed -n "${START_LINE},${END_LINE}p" "${FILE_LIST}")

# Check if we have any files
if [ -z "$JOB_FILES" ]; then
    echo "No files found for job ${JOB_INDEX} (lines ${START_LINE}-${END_LINE})"
    exit 1
fi

# Count actual number of files (filter out empty lines)
NUM_FILES=$(echo "$JOB_FILES" | grep -v '^$' | wc -l)
echo "Job ${JOB_INDEX}: Processing ${NUM_FILES} files (lines ${START_LINE}-${END_LINE})"

# Create temporary file list for Python
TEMP_FILE_LIST="${SCRIPT_DIR}/filelist_job${JOB_INDEX}.txt"
echo "$JOB_FILES" | grep -v '^$' > "${TEMP_FILE_LIST}"

# Create temporary config file
cat > "${TEMP_CONFIG}" << 'CONFIG_EOF'
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

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

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(
CONFIG_EOF

# Append file list to config
while IFS= read -r line; do
    if [ -n "$line" ]; then
        echo "        '${line}'," >> "${TEMP_CONFIG}"
    fi
done < "${TEMP_FILE_LIST}"

# Remove trailing comma from last line and close the list
sed -i '$ s/,$//' "${TEMP_CONFIG}"

# Continue with rest of config
cat >> "${TEMP_CONFIG}" << 'CONFIG_EOF'
    )
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('141X_dataRun3_Express_v3')

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
        massWindow=cms.double(0.025)
    ),
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.025)
    )
)
process.ChiC2To4K = _ChiC2To4K.clone()
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
        massWindow=cms.double(0.090)
    ),
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.090)
    )
)
process.ChiCTo2Ka = _ChiCTo2Ka.clone()
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
        massWindow=cms.double(0.080)
    ),
    cms.PSet(
        name=cms.string('ChiC2'),
        pdgId=cms.int32(445),
        mass=cms.double(3.5562),
        massWindow=cms.double(0.080)
    )
)

process.KshortProducer = _KshortProducer.clone()
process.KshortProducer.trackRecoAlgorithm = cms.InputTag('generalTracks')
process.KshortProducer.vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices')
process.KshortProducer.tkChi2Cut = cms.double(7.0)
process.KshortProducer.tkNhitsCut = cms.int32(3)
process.KshortProducer.tkPtCut = cms.double(0.5)
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
            massWindow=cms.double(0.050)
        )
    )
)
process.ChiC2FromKshortsSequence = cms.Sequence(process.KshortProducer * process.ChiC2FromKshorts)

process.DiKaonWideMass = _DiKaonProducer.clone(
    phiMassCut=cms.double(0.7),
    mKKCutMin=cms.double(0.4),
    mKKCutMax=cms.double(1.5)
)
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
        massWindow=cms.double(0.060)
    )
)
process.ChiC2FromDiKaonPairsSequence = cms.Sequence(process.DiKaonWideMass * process.ChiC2FromDiKaonPairs)

process.ChiCNtuplizer = _ChiCNtuplizer.clone(
    treeName=cms.untracked.string('ChiCNtuple'),
    storeDaughterInfo=cms.untracked.bool(True),
    sources=cms.VPSet(
        cms.PSet(
            name=cms.string('ChiC0_4Pi'),
            pdgId=cms.int32(10441),
            collection=cms.InputTag('ChiCTo4Pi', 'ChiC0')
        ),
        cms.PSet(
            name=cms.string('ChiC2_4Pi'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCTo4Pi', 'ChiC2')
        ),
        cms.PSet(
            name=cms.string('ChiC0_4K'),
            pdgId=cms.int32(10441),
            collection=cms.InputTag('ChiC2To4K', 'ChiC0')
        ),
        cms.PSet(
            name=cms.string('ChiC2_4K'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiC2To4K', 'ChiC2')
        ),
        cms.PSet(
            name=cms.string('ChiC0_2Ka'),
            pdgId=cms.int32(10441),
            collection=cms.InputTag('ChiCTo2Ka', 'ChiC0')
        ),
        cms.PSet(
            name=cms.string('ChiC2_2Ka'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiCTo2Ka', 'ChiC2')
        ),
        cms.PSet(
            name=cms.string('ChiC2_KshortKshort'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiC2FromKshorts', 'ChiC2')
        ),
        cms.PSet(
            name=cms.string('ChiC2_KK'),
            pdgId=cms.int32(445),
            collection=cms.InputTag('ChiC2FromDiKaonPairs', 'ChiC2')
        )
    )
)

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
    fileName=cms.string('${OUTPUT_FILE}'),
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
process.options.numberOfThreads = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.chicNtuple = cms.EndPath(process.ChiCNtuplizer)
process.schedule.append(process.chicNtuple)
CONFIG_EOF

# Clean up temporary file list
rm -f "${TEMP_FILE_LIST}"

# Setup CMSSW environment
cd /afs/cern.ch/work/s/soohwan/private/Analysis/DmesonAna/2025OxygenAnalysis/CMSSW_15_0_9_patch4/src
eval `scramv1 runtime -sh`

# Run the job
echo "Starting job ${JOB_INDEX} at $(date)"
cmsRun "${TEMP_CONFIG}" 2>&1 | tee "${OUTPUT_DIR}/job${JOB_INDEX}.log"

EXIT_CODE=$?

# Clean up temporary config
rm -f "${TEMP_CONFIG}"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Job ${JOB_INDEX} completed successfully at $(date)"
    # Check if output file exists and has reasonable size
    if [ -f "${OUTPUT_FILE}" ] && [ -s "${OUTPUT_FILE}" ]; then
        echo "Output file created: ${OUTPUT_FILE}"
    else
        echo "WARNING: Output file ${OUTPUT_FILE} is missing or empty"
    fi
else
    echo "Job ${JOB_INDEX} failed with exit code ${EXIT_CODE} at $(date)"
fi

exit $EXIT_CODE
