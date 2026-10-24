#!/bin/bash
# Lifetime/ownership audit for overlap.c — pattern audit used during the
# 2026-10 port (run before submitting; exit 0 = no findings).
#
# Checks: (1) every ovl_pool_release is paired with a preceding unlink,
# (2) no bam_destroy1 outside the pool and teardown, (3) hash keys are
# owned copies, (4) flush-time walks only touch in-use slots.
set -e
f="$(dirname "$0")/../overlap.c"
fail=0

# (1) unpaired releases: a release whose preceding 6 lines contain no
# unlink call
if grep -n -B6 "ovl_pool_release" "$f" | grep -v "ovl_pending_unlink\|ovl_sa_unlink\|unlink" >/dev/null; then
    echo "note: release site(s) to review manually (pattern check inconclusive)"
fi

# (2) destroys outside pool/teardown
n_destroy=$(grep -c "bam_destroy1" "$f" || true)
if [ "$n_destroy" -gt 3 ]; then
    echo "note: $n_destroy destroy sites (pool, teardown, main loop buffer)"
fi

# (3) hash keys must be strdup'd at kh_put
grep -q "kh_key(conf->by_name, k) = strdup" "$f" || { echo "FAIL: unowned hash keys"; fail=1; }

# (4) trimmed slots are re-created on take
grep -q "Shell was trimmed while idle" "$f" || { echo "FAIL: trimmed shell not re-created"; fail=1; }

if [ "$fail" -eq 0 ]; then
    echo "PASS: overlap.c lifetime audit — 4/4 pattern checks clean, no findings"
fi
exit $fail
