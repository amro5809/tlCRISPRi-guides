#!/usr/bin/env bash
# run_test.sh — Run design_guides.py on the 3-gene test CSV.
#
# Usage (from repo root or from test/):
#   ./test/run_test.sh /path/to/genome.gbk
#   cd test && ./run_test.sh /path/to/genome.gbk
#
# The first argument must be the path to an annotated GenBank file
# (.gbk / .gbff) for R. palustris CGA009 (or any compatible genome).

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve the repo root regardless of where the script is called from
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

CSV_FILE="${REPO_ROOT}/test/crRNA_test_generation.csv"
SCRIPT="${REPO_ROOT}/design_guides.py"
OUTPUT="${REPO_ROOT}/test/guides_output.csv"

# ---------------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------------
GENBANK_PATH="${1:-}"

if [[ -z "${GENBANK_PATH}" ]]; then
    echo ""
    echo "ERROR: No GenBank file path provided."
    echo ""
    echo "Usage:"
    echo "  ${BASH_SOURCE[0]} /path/to/genome.gbk"
    echo ""
    echo "The genome file must be an annotated GenBank (.gbk / .gbff) file,"
    echo "NOT a plain FASTA (.fna) file, because design_guides.py requires"
    echo "CDS feature annotations to locate RBS and coding regions."
    echo ""
    echo "Example (R. palustris CGA009 from NCBI):"
    echo "  ${BASH_SOURCE[0]} ~/downloads/GCF_000013785.1_genomic.gbff"
    echo ""
    exit 1
fi

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if [[ ! -f "${GENBANK_PATH}" ]]; then
    echo "ERROR: GenBank file not found: ${GENBANK_PATH}"
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: design_guides.py not found at: ${SCRIPT}"
    exit 1
fi

if [[ ! -f "${CSV_FILE}" ]]; then
    echo "ERROR: Test CSV not found at: ${CSV_FILE}"
    exit 1
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "======================================================"
echo "  tlCRISPRi guide design — test run"
echo "======================================================"
echo "  Repo root : ${REPO_ROOT}"
echo "  CSV       : ${CSV_FILE}"
echo "  GenBank   : ${GENBANK_PATH}"
echo "  Output    : ${OUTPUT}"
echo "======================================================"
echo ""

if python "${SCRIPT}" "${CSV_FILE}" "${GENBANK_PATH}" -o "${OUTPUT}"; then
    echo ""
    echo "======================================================"
    echo "  SUCCESS — guides written to:"
    echo "  ${OUTPUT}"
    echo "======================================================"
else
    echo ""
    echo "======================================================"
    echo "  FAILURE — design_guides.py exited with an error."
    echo "  Review the output above for details."
    echo "======================================================"
    exit 1
fi
