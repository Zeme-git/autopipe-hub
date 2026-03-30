#!/bin/bash
# Entrypoint: remove stale Snakemake locks before running
# The pipeline source dir is bind-mounted, so locks from killed runs persist
rm -rf /pipeline/.snakemake/locks 2>/dev/null || true
exec "$@"