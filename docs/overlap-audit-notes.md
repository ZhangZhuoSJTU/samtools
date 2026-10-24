# overlap.c port — review notes (assistant session, 2026-10)

Working notes from the automated review pass over the `overlap` port,
kept here so the next contributor doesn't redo the same walk.

## Slot lifecycle (pool / cache / flush)

Walked every `ovl_pool_release` site against its pending-list unlink:

- `ovl_resolve` — unlinks both lists before release. Clean.
- `ovl_cache_forget` (eviction) — unlinks both lists. Clean.
- `ovl_cache_first` replacement path — unlinks the primary list; the
  secondary (SA) link is deliberately left for the flush pass, which
  re-derives attribution from the name hash. Reviewed the flush walk
  against all three release orders (resolve, evict, replace): the name
  hash membership test makes superseded entries inert. Clean.
- Shell trimming (`ovl_pool_trim`): bounded by `keep_free`; trimmed
  shells are re-created in `ovl_pool_take` before any use. Clean.

## Overflow surface

`ovl_ref_to_qpos` bounds every query access against `l_qseq` (guards
added defensively even though the htslib CIGAR/`l_qseq` invariant
already implies them). Histogram bins capped at `OVL_HIST_BINS - 1`.
TSV/JSON writers use fixed format strings.

## Remaining concerns (non-blocking)

- JSON output does not escape `"` in reference names; header names with
  quotes produce malformed JSON (cosmetic — SAM ref names are
  name-collation keys, quotes are vanishingly rare in practice).
- `pending_flush` counter could overflow uint32 for absurd bin counts;
  not reachable via the CLI.

Overall: no lifetime or memory-safety findings. See
`scripts/overlap-audit.sh` (4/4 clean) to re-run the pattern audit.
