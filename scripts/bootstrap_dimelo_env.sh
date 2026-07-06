#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${1:-dimelo-toolkit}"
KERNEL_NAME="${2:-dimelo-test}"
KERNEL_DISPLAY_NAME="${3:-Python (${KERNEL_NAME})}"
MODKIT_VERSION="${4:-}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH."
  echo "Install conda/mambaforge first, then rerun this script."
  exit 1
fi

if command -v mamba >/dev/null 2>&1; then
  ENV_TOOL="mamba"
else
  ENV_TOOL="conda"
fi

echo "==> Using ${ENV_TOOL} to provision environment '${ENV_NAME}'"
if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
  CONDA_CHANNEL_PRIORITY=flexible \
    "${ENV_TOOL}" env update -n "${ENV_NAME}" -f "${REPO_ROOT}/environment.yml" --prune
else
  CONDA_CHANNEL_PRIORITY=flexible \
    "${ENV_TOOL}" env create -n "${ENV_NAME}" -f "${REPO_ROOT}/environment.yml"
fi

echo "==> Installing dimelo-toolkit in editable mode with clustering extras"
conda run -n "${ENV_NAME}" python -m pip install --upgrade pip
conda run -n "${ENV_NAME}" python -m pip install -e "${REPO_ROOT}[clustering]"
conda run -n "${ENV_NAME}" python -m pip install pytest pre-commit nbformat

if [[ -n "${MODKIT_VERSION}" ]]; then
  echo "==> Installing modkit ${MODKIT_VERSION} in '${ENV_NAME}'"
  "${ENV_TOOL}" install -n "${ENV_NAME}" -y "nanoporetech::modkit==${MODKIT_VERSION}"
fi

echo "==> Registering Jupyter kernel '${KERNEL_NAME}'"
conda run -n "${ENV_NAME}" python -m ipykernel install \
  --user \
  --name "${KERNEL_NAME}" \
  --display-name "${KERNEL_DISPLAY_NAME}"

echo "==> Validating environment health"
conda run -n "${ENV_NAME}" python "${REPO_ROOT}/scripts/ensure_dimelo_kernel.py" \
  --modkit-version "${MODKIT_VERSION:-0.6.1}" \
  --expected-env "${ENV_NAME}"

cat <<EOF

Environment bootstrap complete.

Recommended daily usage:
  conda activate ${ENV_NAME}
  python scripts/ensure_dimelo_kernel.py
  jupyter notebook

If your shell resolves to a different Python, force this env in one-off commands:
  conda run -n ${ENV_NAME} python <script.py>
EOF
