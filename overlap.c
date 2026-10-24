/* overlap.c -- samtools overlap subcommand

    Copyright (C) 2026 Genome Research Ltd.

    Author: yp2838

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/* Reports overlap statistics for read pairs, and optionally repairs
 * inconsistent mate tags (MC, RNEXT/PNEXT mismatches) by emitting a
 * corrected output file.
 *
 * The input does not need to be name-collated: first reads are held in
 * a name-keyed cache of recycled record slots until their mates arrive.
 * The cache is bounded (--max-cache); when full, the least recently
 * inserted entries are evicted and reported as unresolved.
 *
 * Per-reference accounting keeps a secondary list of the pending first
 * reads per reference so that end-of-input reporting can attribute
 * unresolved pairs to the right sequence.
 */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <inttypes.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "samtools.h"
#include "sam_opts.h"

#define DEFAULT_MAX_CACHE 1000000
#define DEFAULT_MIN_OVERLAP 1

KHASH_MAP_INIT_STR(mname, int32_t)

/* One cached first read.  Owned by the slot pool; the name hash maps a
 * QNAME to the pool index of its entry, and the per-reference pending
 * list threads the same indices together for flush-time reporting. */
typedef struct {
    bam1_t *b;          /* the record itself (recycled) */
    uint64_t serial;    /* insertion order, for LRU eviction */
    int32_t next_tid;   /* per-reference pending list link */
    bool in_use;
} ovl_slot_t;

/* Recycled record pool: free slots are threaded on free_next. */
typedef struct {
    ovl_slot_t *slots;
    int32_t n, m;       /* used, allocated */
    int32_t free_head;  /* first free slot, -1 if none */
} ovl_pool_t;

/* Aggregate overlap accounting, per reference and overall. */
typedef struct {
    uint64_t n_pairs;       /* pairs resolved */
    uint64_t n_overlapping; /* pairs with overlap >= min */
    uint64_t sum_overlap;   /* total overlapping bases */
    uint64_t sum_abs_tlen;  /* for consistency reporting */
    uint32_t hist[16];      /* overlap magnitude histogram, log2 bins */
} ovl_stats_t;

typedef struct {
    samFile *in;
    samFile *out;           /* repaired records, or NULL */
    sam_hdr_t *hdr;
    hts_itr_t *iter;
    ovl_pool_t pool;
    khash_t(mname) *by_name;
    int32_t *ref_pending;   /* per-tid head of pending list, -1 if none */
    int n_ref;
    ovl_stats_t *per_ref;   /* sized n_ref */
    ovl_stats_t total;
    uint64_t serial;
    int64_t max_cache;
    int min_overlap;
    int min_mapq;
    int flag_on, flag_off;
    bool tsv;
    bool per_ref_report;
    bool repair;
    int64_t n_evicted, n_unresolved, n_singletons, n_records;
} ovl_conf_t;

static void ovl_stats_add(ovl_stats_t *st, int64_t overlap, int64_t tlen,
                          int min_overlap) {
    st->n_pairs++;
    if (overlap >= min_overlap && overlap >= 1) {
        st->n_overlapping++;
        st->sum_overlap += overlap;
        int bin = 0;
        while (bin < 15 && (overlap >>= 1) > 0)
            bin++;
        st->hist[bin]++;
    }
    st->sum_abs_tlen += tlen < 0 ? -tlen : tlen;
}

/* Return the number of reference bases shared by two mapped records on
 * the same reference.  Both intervals are half-open [pos, end). */
static int64_t ovl_pair_overlap(const bam1_t *a, const bam1_t *b) {
    hts_pos_t a_start = a->core.pos, a_end = bam_endpos(a);
    hts_pos_t b_start = b->core.pos, b_end = bam_endpos(b);
    hts_pos_t lo = a_start > b_start ? a_start : b_start;
    hts_pos_t hi = a_end < b_end ? a_end : b_end;
    return hi > lo ? hi - lo : 0;
}

/* Mate cigar as a string, for MC tag maintenance.  Writes into ks. */
static int ovl_cigar_kstring(const bam1_t *b, kstring_t *ks) {
    uint32_t *cigar = bam_get_cigar(b);
    ks->l = 0;
    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        if (ksprintf(ks, "%"PRIu32"%c", bam_cigar_oplen(cigar[i]),
                     bam_cigar_opchr(cigar[i])) < 0)
            return -1;
    }
    return 0;
}

/* True when a record's mate coordinates disagree with the mate itself. */
static bool ovl_mate_tags_stale(const bam1_t *r1, const bam1_t *r2) {
    if (r1->core.mtid != r2->core.tid) return true;
    if (r1->core.mpos != r2->core.pos) return true;
    /* TLEN sign must live on the leftmost end; a mismatch is stale. */
    int64_t span = bam_endpos(r2) - r1->core.pos;
    if (r1->core.isize > 0) return r1->core.isize != span;
    if (r2->core.isize < 0) return r2->core.isize != -span;
    return true;
}

/* Rewrite r1's mate bookkeeping from the (newly seen) r2, and vice
 * versa.  Also refreshes the MC tags when present. */
static int ovl_repair_pair(bam1_t *r1, bam1_t *r2) {
    kstring_t ks = KS_INITIALIZE;
    int rc = 0;

    r1->core.mtid = r2->core.tid;
    r1->core.mpos = r2->core.pos;
    r2->core.mtid = r1->core.tid;
    r2->core.mpos = r1->core.pos;

    int64_t span = bam_endpos(r2) - r1->core.pos;
    if (span > 0) {
        r1->core.isize = span;
        r2->core.isize = -span;
    } else {
        r1->core.isize = bam_endpos(r1) - r2->core.pos;
        r2->core.isize = -r1->core.isize;
    }

    if (bam_aux_get(r1, "MC") || bam_aux_get(r2, "MC")) {
        if (ovl_cigar_kstring(r2, &ks) < 0) { rc = -1; goto out; }
        if (bam_aux_update_str(r1, "MC", ks.l, ks.s) < 0) { rc = -1; goto out; }
        if (ovl_cigar_kstring(r1, &ks) < 0) { rc = -1; goto out; }
        if (bam_aux_update_str(r2, "MC", ks.l, ks.s) < 0) { rc = -1; goto out; }
    }
 out:
    ks_free(&ks);
    return rc;
}

/* --- slot pool ------------------------------------------------- */

static int ovl_pool_init(ovl_pool_t *p) {
    p->slots = NULL;
    p->n = p->m = 0;
    p->free_head = -1;
    return 0;
}

/* Take a slot from the free list, or append a fresh one.  Returns the
 * slot index, or -1 on allocation failure. */
static int32_t ovl_pool_take(ovl_pool_t *p) {
    if (p->free_head >= 0) {
        int32_t s = p->free_head;
        p->free_head = p->slots[s].next_tid; /* free list reuses the link */
        p->slots[s].in_use = true;
        return s;
    }
    if (p->n == p->m) {
        int32_t m = p->m ? p->m * 2 : 1024;
        ovl_slot_t *ns = realloc(p->slots, m * sizeof(*ns));
        if (!ns) return -1;
        p->slots = ns;
        p->m = m;
    }
    int32_t s = p->n++;
    p->slots[s].b = bam_init1();
    if (!p->slots[s].b) return -1;
    p->slots[s].in_use = true;
    p->slots[s].next_tid = -1;
    return s;
}

/* Hand a slot back to the free list; the record itself is kept for
 * reuse.  The caller must already have unlinked it from the pending
 * list (the link field carries the free-list chain from here on). */
static void ovl_pool_release(ovl_pool_t *p, int32_t s) {
    p->slots[s].in_use = false;
    p->slots[s].next_tid = p->free_head;
    p->free_head = s;
}

static void ovl_pool_destroy(ovl_pool_t *p) {
    for (int32_t i = 0; i < p->n; i++)
        if (p->slots[i].b)
            bam_destroy1(p->slots[i].b);
    free(p->slots);
    p->slots = NULL;
    p->n = p->m = 0;
    p->free_head = -1;
}

/* --- cache operations ------------------------------------------ */

/* Remove slot s from the per-reference pending list it is on.  The
 * pending lists are singly linked in insertion order, so unlinking
 * walks from the head; pending lists are short in practice because
 * coordinate-sorted input resolves pairs within a small window. */
static void ovl_pending_unlink(ovl_conf_t *conf, int32_t s) {
    int tid = conf->pool.slots[s].b->core.tid;
    if (tid < 0 || tid >= conf->n_ref)
        return;
    int32_t *cur = &conf->ref_pending[tid];
    while (*cur >= 0) {
        if (*cur == s) {
            *cur = conf->pool.slots[s].next_tid;
            conf->pool.slots[s].next_tid = -1;
            return;
        }
        cur = &conf->pool.slots[*cur].next_tid;
    }
}

static void ovl_pending_push(ovl_conf_t *conf, int32_t s) {
    int tid = conf->pool.slots[s].b->core.tid;
    if (tid < 0 || tid >= conf->n_ref)
        return;
    int32_t *cur = &conf->ref_pending[tid];
    while (*cur >= 0)
        cur = &conf->pool.slots[*cur].next_tid;
    *cur = s;
}

/* Drop a cached first read: unlink from pending, forget the name, and
 * recycle the slot. */
static void ovl_cache_forget(ovl_conf_t *conf, int32_t s, const char *name) {
    ovl_pending_unlink(conf, s);
    khint_t k = kh_get(mname, conf->by_name, name);
    if (k != kh_end(conf->by_name)) {
        free((char *) kh_key(conf->by_name, k));
        kh_del(mname, conf->by_name, k);
    }
    ovl_pool_release(&conf->pool, s);
    conf->n_evicted++;
}

/* Resolve a completed pair: stats, optional repair and output. */
static int ovl_emit_pair(ovl_conf_t *conf, bam1_t *r1, bam1_t *r2) {
    int tid = r1->core.tid;
    int64_t ov = 0;
    bool same_ref = tid >= 0 && tid == r2->core.tid;
    if (same_ref)
        ov = ovl_pair_overlap(r1, r2);
    {
        ovl_stats_t *st = tid >= 0 && tid < conf->n_ref
            ? &conf->per_ref[tid] : NULL;
        if (st) ovl_stats_add(st, ov, r1->core.isize, conf->min_overlap);
        ovl_stats_add(&conf->total, ov, r1->core.isize, conf->min_overlap);
    }
    if (conf->repair && same_ref) {
        if (ovl_mate_tags_stale(r1, r2)) {
            if (ovl_repair_pair(r1, r2) < 0)
                return -1;
        }
        if (sam_write1(conf->out, conf->hdr, r1) < 0) return -1;
        if (sam_write1(conf->out, conf->hdr, r2) < 0) return -1;
    }
    return 0;
}

/* A second read arrived: pair it with its cached first read if we
 * still hold it.  Returns 0 on success, -1 on error. */
static int ovl_resolve(ovl_conf_t *conf, bam1_t *r2) {
    const char *name = bam_get_qname(r2);
    khint_t k = kh_get(mname, conf->by_name, name);
    if (k == kh_end(conf->by_name)) {
        conf->n_singletons++;
        return 0;
    }
    int32_t s = kh_value(conf->by_name, k);
    bam1_t *r1 = conf->pool.slots[s].b;

    if (ovl_emit_pair(conf, r1, r2) < 0)
        return -1;

    /* Resolve in place: forget the name and recycle the slot, but keep
     * the pair attribution on the pending list until flush; the unlink
     * must happen before release, as release reuses the link field. */
    ovl_pending_unlink(conf, s);
    free((char *) kh_key(conf->by_name, k));
    kh_del(mname, conf->by_name, k);
    ovl_pool_release(&conf->pool, s);
    return 0;
}

/* Cache a first read.  Under cache pressure, evict the oldest entries
 * until there is room.  A duplicate name replaces the older copy. */
static int ovl_cache_first(ovl_conf_t *conf, bam1_t *b) {
    const char *name = bam_get_qname(b);
    khint_t k = kh_get(mname, conf->by_name, name);

    if (k != kh_end(conf->by_name)) {
        /* Same template name seen again before its mate: the previous
         * copy is superseded (e.g. supplementary alignment). */
        int32_t old = kh_value(conf->by_name, k);
        ovl_pending_unlink(conf, old);
        ovl_pool_release(&conf->pool, old);
        k = kh_end(conf->by_name);
    }

    while ((int64_t) kh_size(conf->by_name) >= conf->max_cache) {
        /* Evict the least recently inserted entry. */
        uint64_t best = UINT64_MAX;
        int32_t victim = -1;
        for (k = kh_begin(conf->by_name); k != kh_end(conf->by_name); k++) {
            if (!kh_exist(conf->by_name, k)) continue;
            int32_t s = kh_value(conf->by_name, k);
            if (conf->pool.slots[s].serial < best) {
                best = conf->pool.slots[s].serial;
                victim = s;
            }
        }
        if (victim < 0)
            break;
        ovl_cache_forget(conf, victim, bam_get_qname(conf->pool.slots[victim].b));
    }

    int32_t s = ovl_pool_take(&conf->pool);
    if (s < 0)
        return -1;
    if (bam_copy1(conf->pool.slots[s].b, b) < 0)
        return -1;
    conf->pool.slots[s].serial = ++conf->serial;

    int ret;
    k = kh_put(mname, conf->by_name, name, &ret);
    if (ret < 0)
        return -1;
    if (ret > 0)
        kh_key(conf->by_name, k) = strdup(name);
    if (ret == 0) {
        /* Name reappeared while we held it: keep the newer copy. */
        int32_t prev = kh_value(conf->by_name, k);
        ovl_pending_unlink(conf, prev);
        ovl_pool_release(&conf->pool, prev);
    }
    kh_value(conf->by_name, k) = s;

    ovl_pending_push(conf, s);
    return 0;
}

/* --- reporting -------------------------------------------------- */

static void ovl_hist_row(FILE *fp, const ovl_stats_t *st) {
    uint64_t total_h = 0;
    for (int i = 0; i < 16; i++) total_h += st->hist[i];
    for (int i = 0; i < 16; i++) {
        if (st->hist[i] == 0) continue;
        fprintf(fp, "\t2^%02d-%02d\t%u\t%g%%\n", i, i + 1, st->hist[i],
                100.0 * st->hist[i] / (total_h ? total_h : 1));
    }
}

static void ovl_hash_destroy(ovl_conf_t *conf) {
    if (!conf->by_name)
        return;
    khint_t k;
    for (k = kh_begin(conf->by_name); k != kh_end(conf->by_name); k++)
        if (kh_exist(conf->by_name, k))
            free((char *) kh_key(conf->by_name, k));
    kh_destroy(mname, conf->by_name);
    conf->by_name = NULL;
}

static void ovl_report(FILE *fp, const ovl_conf_t *conf) {
    const ovl_stats_t *t = &conf->total;
    if (conf->tsv) {
        fprintf(fp, "#rname\tpairs\toverlapping\tsum_overlap\tmean_overlap\n");
        if (conf->per_ref_report) {
            for (int i = 0; i < conf->n_ref; i++) {
                if (conf->per_ref[i].n_pairs == 0) continue;
                fprintf(fp, "%s\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%g\n",
                        sam_hdr_tid2name(conf->hdr, i),
                        conf->per_ref[i].n_pairs,
                        conf->per_ref[i].n_overlapping,
                        conf->per_ref[i].sum_overlap,
                        conf->per_ref[i].n_overlapping
                            ? conf->per_ref[i].sum_overlap
                              / (double) conf->per_ref[i].n_overlapping : 0);
            }
        }
        fprintf(fp, "*\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%g\n",
                t->n_pairs, t->n_overlapping, t->sum_overlap,
                t->n_overlapping ? t->sum_overlap / (double) t->n_overlapping : 0);
    } else {
        fprintf(fp, "#overlapping pairs\t%"PRIu64"/%"PRIu64"\n",
                t->n_overlapping, t->n_pairs);
        fprintf(fp, "#overlapping bases\t%"PRIu64"\n", t->sum_overlap);
        fprintf(fp, "#mean overlap\t%g\n",
                t->n_overlapping ? t->sum_overlap / (double) t->n_overlapping : 0);
        fprintf(fp, "#mean |TLEN|\t%g\n",
                t->n_pairs ? t->sum_abs_tlen / (double) t->n_pairs : 0);
    }
    if (conf->per_ref_report && !conf->tsv) {
        fprintf(fp, "#overlap histogram (all references)\n");
        ovl_hist_row(fp, t);
    }
    fprintf(fp, "#records\t%"PRId64"\n", conf->n_records);
    fprintf(fp, "#unresolved firsts\t%"PRId64"\n", conf->n_unresolved);
    fprintf(fp, "#mateless seconds\t%"PRId64"\n", conf->n_singletons);
    fprintf(fp, "#evictions\t%"PRId64"\n", conf->n_evicted);
}

/* --- option parsing and main loop ------------------------------- */

static int usage(FILE *fp) {
    fprintf(fp,
"Usage: samtools overlap [options] in.sam|in.bam|in.cram\n"
"\n"
"Options:\n"
"  -o FILE       write pairs with repaired mate tags to FILE\n"
"  -r REG        restrict to region (requires an index)\n"
"  -c INT        minimum overlapping bases for a pair to count [1]\n"
"  -q INT        minimum mapping quality [0]\n"
"  -m, --max-cache INT   cached first reads before eviction [%d]\n"
"  -t, --tsv     TSV output instead of the default summary\n"
"  -P, --per-ref report per reference (and a histogram)\n"
"  -f, --require-flags FLAG  skip reads without any of FLAG\n"
"  -F, --excl-flags FLAG    skip reads with any of FLAG\n"
"                [UNMAP,SECONDARY,QCFAIL,DUP]\n",
            DEFAULT_MAX_CACHE);
    sam_global_opt_help(fp, "-.--.--.");
    return EXIT_SUCCESS;
}

int main_overlap(int argc, char *argv[]) {
    int ret = EXIT_FAILURE;
    ovl_conf_t conf;
    memset(&conf, 0, sizeof(conf));
    conf.max_cache = DEFAULT_MAX_CACHE;
    conf.min_overlap = DEFAULT_MIN_OVERLAP;
    conf.flag_off = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    FILE *report = stdout;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    const char *out_fn = NULL, *reg = NULL;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '@'),
        {"tsv", no_argument, NULL, 't'},
        {"per-ref", no_argument, NULL, 'P'},
        {"repair", no_argument, NULL, 'R'},
        {"max-cache", required_argument, NULL, 'm'},
        {"min-overlap", required_argument, NULL, 'c'},
        {"min-mapq", required_argument, NULL, 'q'},
        {"require-flags", required_argument, NULL, 'f'},
        {"excl-flags", required_argument, NULL, 'F'},
        {"exclude-flags", required_argument, NULL, 'F'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:r:c:q:f:F:m:tP", lopts, NULL)) >= 0) {
        switch (c) {
        case 'o': out_fn = optarg; break;
        case 'r': reg = optarg; break;
        case 'c': conf.min_overlap = atoi(optarg); break;
        case 'q': conf.min_mapq = atoi(optarg); break;
        case 'm': conf.max_cache = strtoll(optarg, NULL, 0); break;
        case 't': conf.tsv = true; break;
        case 'P': conf.per_ref_report = true; break;
        case 'f':
            if ((conf.flag_on = bam_str2flag(optarg)) < 0) {
                print_error("overlap", "unknown flag \"%s\"", optarg);
                return usage(stderr);
            }
            break;
        case 'F':
            if ((conf.flag_off = bam_str2flag(optarg)) < 0) {
                print_error("overlap", "unknown flag \"%s\"", optarg);
                return usage(stderr);
            }
            break;
        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': return usage(stderr);
        }
    }
    if (argc != optind + 1 || conf.max_cache <= 0)
        return usage(stderr);

    conf.in = sam_open_format(argv[optind], "r", &ga.in);
    if (!conf.in) {
        print_error_errno("overlap", "cannot open \"%s\"", argv[optind]);
        goto out;
    }
    conf.hdr = sam_hdr_read(conf.in);
    if (!conf.hdr) {
        print_error_errno("overlap", "cannot read header from \"%s\"", argv[optind]);
        goto out;
    }
    if (out_fn) {
        conf.repair = true;
        conf.out = sam_open_format(out_fn, "wm", &ga.out);
        if (!conf.out) {
            print_error_errno("overlap", "cannot open \"%s\" for writing", out_fn);
            goto out;
        }
        if (sam_hdr_write(conf.out, conf.hdr) < 0) goto out;
    }
    if (reg) {
        hts_idx_t *idx = sam_index_load(conf.in, argv[optind]);
        if (!idx) {
            print_error_errno("overlap", "cannot load index for \"%s\"", argv[optind]);
            goto out;
        }
        conf.iter = sam_itr_querys(idx, conf.hdr, reg);
        hts_idx_destroy(idx);
        if (!conf.iter) {
            print_error("overlap", "cannot parse region \"%s\"", reg);
            goto out;
        }
    }

    conf.n_ref = sam_hdr_nref(conf.hdr);
    conf.per_ref = calloc(conf.n_ref ? conf.n_ref : 1, sizeof(*conf.per_ref));
    conf.ref_pending = malloc((conf.n_ref ? conf.n_ref : 1) * sizeof(int32_t));
    conf.by_name = kh_init(mname);
    if (!conf.per_ref || !conf.ref_pending || !conf.by_name) {
        print_error_errno("overlap", "out of memory");
        goto out;
    }
    for (int i = 0; i < conf.n_ref; i++)
        conf.ref_pending[i] = -1;
    ovl_pool_init(&conf.pool);

    bam1_t *b = bam_init1();
    if (!b) goto out;
    int r;
    while ((r = conf.iter
                ? sam_itr_next(conf.in, conf.iter, b)
                : sam_read1(conf.in, conf.hdr, b)) >= 0) {
        conf.n_records++;
        if (b->core.flag & conf.flag_off) continue;
        if (conf.flag_on && !(b->core.flag & conf.flag_on)) continue;
        if (b->core.qual < conf.min_mapq) continue;

        bool read1 = (b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREAD2);
        bool read2 = (b->core.flag & BAM_FREAD2) && !(b->core.flag & BAM_FREAD1);
        if (read1) {
            if (ovl_cache_first(&conf, b) < 0) {
                print_error("overlap", "cannot cache record");
                goto out;
            }
        } else if (read2) {
            if (ovl_resolve(&conf, b) < 0) {
                print_error_errno("overlap", "cannot resolve pair");
                goto out;
            }
        } else {
            /* Unpaired or malformed flags: nothing to reconcile. */
            continue;
        }
    }
    if (r < -1) {
        print_error("overlap", "error reading input");
        goto out;
    }
    bam_destroy1(b);
    b = NULL;

    /* End of input: every pending first read is unresolved. */
    for (int i = 0; i < conf.n_ref; i++) {
        int32_t s = conf.ref_pending[i];
        while (s >= 0) {
            conf.n_unresolved++;
            s = conf.pool.slots[s].next_tid;
        }
    }
    ovl_report(report, &conf);
    ret = EXIT_SUCCESS;

 out:
    if (b) bam_destroy1(b);
    if (conf.iter) hts_itr_destroy(conf.iter);
    ovl_hash_destroy(&conf);
    ovl_pool_destroy(&conf.pool);
    free(conf.per_ref);
    free(conf.ref_pending);
    if (conf.out) sam_close(conf.out);
    if (conf.in) sam_close(conf.in);
    if (conf.hdr) sam_hdr_destroy(conf.hdr);
    sam_global_args_free(&ga);
    return ret;
}
