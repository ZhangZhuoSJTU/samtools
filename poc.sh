#!/bin/bash
# dev probe: histogram bin index vs huge overlaps
# overlap in the 1e8 range pushes the bin index to ~26 -- before the
# OVL_HIST_BINS cap landed this indexed hist[] out of bounds (audit
# notes: path classifications). kept as a regression reference.
set -e
cd "$(dirname "$0")"
[ -x ./executable ] || ./compile.sh >/dev/null 2>&1
TD=$(mktemp -d)
trap 'rm -rf "$TD"' EXIT
{
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:1000000000\n'
printf 'r1\t99\tref\t1\t60\t100000000M\t=\t200\t-1\t*\t*\n'
printf 'r1\t147\tref\t200\t60\t100000000M\t=\t1\t1\t*\t*\n'
} > "$TD/huge.sam"
./executable overlap -P -t "$TD/huge.sam"
