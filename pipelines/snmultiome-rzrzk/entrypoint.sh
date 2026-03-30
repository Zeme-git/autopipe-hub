#!/bin/bash
# Entrypoint: unlock stale Snakemake locks from bind-mounted pipeline dir
# Uses snakemake's built-in unlock mechanism (safe, no destructive commands)
cd /pipeline
conda run --no-capture-output -n pipeline snakemake --unlock 2>/dev/null || true
exec "$@"