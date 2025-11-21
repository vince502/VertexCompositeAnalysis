#!/bin/bash

# Helper script to submit Condor jobs
# Calculates the number of jobs needed based on file count

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILE_LIST="${SCRIPT_DIR}/file.txt"
FILES_PER_JOB=20
OUTPUT_DIR="/eos/cms/store/group/phys_heavyions_ops/soohwan/Run3_2025/Oxy/ChiC"

# Count total files
TOTAL_FILES=$(wc -l < "${FILE_LIST}")
NUM_JOBS=$(( (TOTAL_FILES + FILES_PER_JOB - 1) / FILES_PER_JOB ))

echo "Total files: ${TOTAL_FILES}"
echo "Files per job: ${FILES_PER_JOB}"
echo "Number of jobs needed: ${NUM_JOBS}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Update condor.sub with correct number of jobs and output directory
sed -i "s/^queue.*/queue ${NUM_JOBS}/" "${SCRIPT_DIR}/condor.sub"
# Update output directory in condor.sub arguments line
sed -i "s|arguments = \$(Process) [0-9]* .*|arguments = \$(Process) ${FILES_PER_JOB} ${OUTPUT_DIR}|" "${SCRIPT_DIR}/condor.sub"
# Update output/error/log paths in condor.sub
sed -i "s|output = .*|output = ${OUTPUT_DIR}/job\$(Process).out|" "${SCRIPT_DIR}/condor.sub"
sed -i "s|error = .*|error = ${OUTPUT_DIR}/job\$(Process).err|" "${SCRIPT_DIR}/condor.sub"
sed -i "s|log = .*|log = ${OUTPUT_DIR}/condor.log|" "${SCRIPT_DIR}/condor.sub"

echo ""
echo "Updated condor.sub with:"
echo "  - Number of jobs: ${NUM_JOBS}"
echo "  - Files per job: ${FILES_PER_JOB}"
echo "  - Output directory: ${OUTPUT_DIR}"
echo ""
echo "To submit jobs, run:"
echo "  condor_submit ${SCRIPT_DIR}/condor.sub"
echo ""
echo "Note: Using runCondor.sh which processes files in batch mode"
echo "      Each job processes ${FILES_PER_JOB} files using inputFileList argument"
echo ""
echo "To check job status:"
echo "  condor_q"
echo ""
echo "To check specific job output:"
echo "  tail -f ${OUTPUT_DIR}/job0.out"
echo ""
