#!/bin/bash -ex

# Script to run CMSSW job with a subset of files from file.txt
# Usage: ./runCondor.sh <job_index> <files_per_job> <output_dir>

JOB_INDEX=$1
FILES_PER_JOB=${2:-20}
OUTPUT_DIR=${3:-/eos/cms/store/group/phys_heavyions_ops/soohwan/Run3_2025/Oxy/ChiC}

# Calculate line range for this job
START_LINE=$((JOB_INDEX * FILES_PER_JOB + 1))
END_LINE=$((START_LINE + FILES_PER_JOB - 1))

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INITIAL_DIR="$(pwd)"

# Config file: check current directory first (Condor transfers files here), then SCRIPT_DIR
# Use absolute paths to avoid issues when changing directories later
if [ -f "${INITIAL_DIR}/pp2024_ChiC_MB_cfg_v1_ForCondor.py" ]; then
    CONFIG_FILE="${INITIAL_DIR}/pp2024_ChiC_MB_cfg_v1_ForCondor.py"
elif [ -f "./pp2024_ChiC_MB_cfg_v1_ForCondor.py" ]; then
    CONFIG_FILE="$(cd "$(dirname "./pp2024_ChiC_MB_cfg_v1_ForCondor.py")" && pwd)/$(basename "./pp2024_ChiC_MB_cfg_v1_ForCondor.py")"
elif [ -f "${SCRIPT_DIR}/pp2024_ChiC_MB_cfg_v1_ForCondor.py" ]; then
    CONFIG_FILE="${SCRIPT_DIR}/pp2024_ChiC_MB_cfg_v1_ForCondor.py"
else
    echo "ERROR: Config file pp2024_ChiC_MB_cfg_v1_ForCondor.py not found!"
    echo "  Checked: ${INITIAL_DIR}/pp2024_ChiC_MB_cfg_v1_ForCondor.py"
    echo "  Checked: ${SCRIPT_DIR}/pp2024_ChiC_MB_cfg_v1_ForCondor.py"
    exit 1
fi

# File list: check current directory first (Condor transfers files here), then SCRIPT_DIR
if [ -f "${INITIAL_DIR}/file.txt" ]; then
    FILE_LIST="${INITIAL_DIR}/file.txt"
elif [ -f "./file.txt" ]; then
    FILE_LIST="$(cd "$(dirname "./file.txt")" && pwd)/$(basename "./file.txt")"
elif [ -f "${SCRIPT_DIR}/file.txt" ]; then
    FILE_LIST="${SCRIPT_DIR}/file.txt"
else
    echo "ERROR: File list file.txt not found!"
    echo "  Checked: ${INITIAL_DIR}/file.txt"
    echo "  Checked: ${SCRIPT_DIR}/file.txt"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

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

# Create temporary file list for batch processing
TEMP_FILE_LIST="${OUTPUT_DIR}/filelist_job${JOB_INDEX}.txt"
echo "$JOB_FILES" | grep -v '^$' > "${TEMP_FILE_LIST}"

# Output file
OUTPUT_FILE="${OUTPUT_DIR}/chic_combination_tree_job${JOB_INDEX}.root"

# Storage site (can be overridden via environment variable or argument)
STORAGE_SITE=${STORAGE_SITE:-infn}  # Default: cern, options: cern, fnal, local, or custom redirector

# Setup CMSSW environment
CMSSW_BASE="/afs/cern.ch/work/s/soohwan/private/Analysis/DmesonAna/2025OxygenAnalysis/CMSSW_15_0_9_patch4"
cd "${CMSSW_BASE}/src"
eval `scramv1 runtime -sh`

# Run the job with inputFileList argument
echo "Starting job ${JOB_INDEX} at $(date)"
echo "Config file: ${CONFIG_FILE}"
echo "Input file list: ${TEMP_FILE_LIST}"
echo "Output file: ${OUTPUT_FILE}"
echo "Storage site: ${STORAGE_SITE}"
echo "Working directory: $(pwd)"
echo "CMSSW_BASE: ${CMSSW_BASE}"

# Verify config file exists (use absolute path)
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Config file not found at ${CONFIG_FILE}"
    echo "Current directory: $(pwd)"
    echo "Looking for config file..."
    ls -la "${CONFIG_FILE}" 2>&1 || echo "File does not exist"
    exit 1
fi

# Verify temp file list exists
if [ ! -f "${TEMP_FILE_LIST}" ]; then
    echo "ERROR: Temporary file list not found at ${TEMP_FILE_LIST}"
    exit 1
fi

# Run cmsRun with proper arguments
# Use absolute path for config file to avoid any path issues
cmsRun "${CONFIG_FILE}" inputFileList="${TEMP_FILE_LIST}" outputFile="${OUTPUT_FILE}" storageSite="${STORAGE_SITE}" 2>&1 | tee "${OUTPUT_DIR}/job${JOB_INDEX}.log"

EXIT_CODE=$?

# Clean up temporary file list
rm -f "${TEMP_FILE_LIST}"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Job ${JOB_INDEX} completed successfully at $(date)"
    # Check if output file exists and has reasonable size
    if [ -f "${OUTPUT_FILE}" ] && [ -s "${OUTPUT_FILE}" ]; then
        echo "Output file created: ${OUTPUT_FILE}"
        ls -lh "${OUTPUT_FILE}"
    else
        echo "WARNING: Output file ${OUTPUT_FILE} is missing or empty"
    fi
else
    echo "Job ${JOB_INDEX} failed with exit code ${EXIT_CODE} at $(date)"
fi

exit $EXIT_CODE
