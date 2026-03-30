#!/bin/bash
# Entrypoint: clear stale Snakemake lock files from bind-mounted pipeline dir
# Only removes the specific lock files, not the entire .snakemake directory
if [ -d /pipeline/.snakemake/locks ]; then
  find /pipeline/.snakemake/locks -type f -delete 2>/dev/null || true
fi
exec "$@"