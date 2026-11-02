#!/bin/bash
# Basic functional checks for `samtools overlap`.
set -e
EXE="${1:-./executable}"
TD=$(mktemp -d)
trap 'rm -rf "$TD"' EXIT
cat > "$TD/in.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:1000
p1	99	ref	10	60	10M	=	15	15	ACGTACGTAC	FFFFFFFFFF
p1	147	ref	15	60	10M	=	10	-15	ACGTACGTAC	FFFFFFFFFF
p2	99	ref	100	60	10M	=	200	110	ACGTACGTAC	FFFFFFFFFF
p2	147	ref	200	60	10M	=	100	-110	ACGTACGTAC	FFFFFFFFFF
SAM
# overlapping pair found, mean overlap 5
out=$("$EXE" overlap "$TD/in.sam")
echo "$out" | grep -q "#overlapping pairs	1/2"
echo "$out" | grep -q "#overlapping bases	5"
echo "$out" | grep -q "#mean overlap	5"
# TSV mode emits per-reference rows with -P
"$EXE" overlap -t -P "$TD/in.sam" | grep -q "^ref	2	1	5	5"
# JSON mode is valid enough to carry the totals
"$EXE" overlap -J "$TD/in.sam" | grep -q '"overlapping": 1'
# identity filter excludes the shifted pair
"$EXE" overlap -I 80 "$TD/in.sam" | grep -q "#overlapping pairs	0/2"
# shuffled input still resolves via the cache
grep -v '^p1	147' "$TD/in.sam" > "$TD/shuf.sam"
echo 'p1	147	ref	15	60	10M	=	10	-15	ACGTACGTAC	FFFFFFFFFF' >> "$TD/shuf.sam"
"$EXE" overlap "$TD/shuf.sam" | grep -q "#overlapping pairs	1/2"
# repair mode writes both mates
"$EXE" overlap -o "$TD/fix.bam" "$TD/in.sam" >/dev/null
"$EXE" view "$TD/fix.bam" | wc -l | grep -q 4
echo "overlap: all checks passed"
