#!/bin/bash
# ProgramBench-style build: must produce ./executable at the repo root.
# samtools needs a version-matched HTSlib source tree; fetch the matching
# release and build it in-tree (the course build has network; the isolated
# eval build does not, but that only runs on the gold image).
set -e
cd "$(dirname "$0")"
HV=1.23.1
if [ ! -d "htslib-$HV" ]; then
    curl -fsSL "https://github.com/samtools/htslib/releases/download/$HV/htslib-$HV.tar.bz2" | tar xj
fi
autoreconf -i
./configure --with-htslib="htslib-$HV" --without-curses
make -j2
cp samtools executable
