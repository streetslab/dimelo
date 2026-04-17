#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-${CONDA_DEFAULT_ENV:-dimelo-toolkit}}"
MODKIT_VERSION="${2:-}"

if [[ -z "${MODKIT_VERSION}" ]]; then
  echo "Usage: bash scripts/update_modkit_version.sh <conda_env_name> <modkit_version>"
  echo "Example: bash scripts/update_modkit_version.sh dimelo-toolkit 0.6.1"
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH."
  exit 1
fi

if command -v mamba >/dev/null 2>&1; then
  ENV_TOOL="mamba"
else
  ENV_TOOL="conda"
fi

echo "==> Installing nanoporetech::modkit==${MODKIT_VERSION} in env '${ENV_NAME}'"
"${ENV_TOOL}" install -n "${ENV_NAME}" -y "nanoporetech::modkit==${MODKIT_VERSION}"

echo "==> Verifying in '${ENV_NAME}'"
conda run -n "${ENV_NAME}" modkit --version
conda run -n "${ENV_NAME}" python scripts/ensure_dimelo_kernel.py \
  --modkit-version "${MODKIT_VERSION}" \
  --expected-env "${ENV_NAME}"

