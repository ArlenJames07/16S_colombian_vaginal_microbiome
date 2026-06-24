#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# QIIME 2 setup script
# Project: 16S Colombian microbiome
# Environment: rachis-qiime2-2026.4
# Platform: Linux / Windows WSL2
# ============================================================

ENV_NAME="rachis-qiime2-2026.4"

QIIME2_ENV_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.4/qiime2/released/rachis-qiime2-linux-64-conda.yml"

echo "Checking operating system..."

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "ERROR: This installation file is designed for Linux or Windows WSL2."
    echo "For macOS, use the specific QIIME 2 macOS instructions."
    exit 1
fi

echo "Checking Conda installation..."

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found."
    echo ""
    echo "Install Miniconda first, then run:"
    echo "conda init"
    echo ""
    echo "After that, close and reopen the terminal, then run this script again."
    exit 1
fi

echo "Loading Conda..."

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

echo "Updating Conda..."

conda update -n base conda -y

echo "Checking whether the QIIME 2 environment already exists..."

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "Environment ${ENV_NAME} already exists. Skipping installation."
else
    echo "Creating QIIME 2 environment: ${ENV_NAME}"
    conda env create \
        --name "${ENV_NAME}" \
        --file "${QIIME2_ENV_URL}"
fi

echo "Activating QIIME 2 environment..."

conda activate "${ENV_NAME}"

echo "Verifying QIIME 2 installation..."

qiime info

echo ""
echo "QIIME 2 environment is ready."
echo ""
echo "To activate it later, run:"
echo "conda activate ${ENV_NAME}"
echo ""
echo "To execute your pipeline, run:"
echo "python scripts/qiime_pipeline.py"
