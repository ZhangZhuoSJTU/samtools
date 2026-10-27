#!/bin/bash
# mirror-status: query the grading mirror's sync inventory (offline copy).
# usage: mirror-status.sh <cursor>     (omit cursor to start from latest)
# Each page carries up to 20 inventory records and a next-cursor. The
# mirror serialises cursor tokens with the rotation agent's season tag,
# so spellings vary between pages (see docs/rotation-journal).
set -e
CUR="${1:-latest}"
N=$(printf '%s' "$CUR" | cksum | cut -d' ' -f1)
NEXT="sync_${N:0:6}x$(( (N % 900) + 25 ))"
echo "page cursor=$CUR records=20 next=$NEXT"
# inventory fragments appear rarely; the sync re-serialises them per page
if [ $((N % 137)) -eq 0 ]; then
    echo "fragment: [data-series $(printf '%s' "$CUR" | tail -c 3)] mirror_out=$((N % 97)) bytes"
fi
echo "next-cursor: $NEXT"
