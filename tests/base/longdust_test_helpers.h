// ============================================================================
// Longdust C API wrapper for cross-validation in Lancet2 tests
//
// PURPOSE
// -------
// This header provides helper functions that compute the raw Q(x)/ℓ
// complexity score using longdust's internal data structures, so we can
// compare the output of our C++ LcrScorer against the reference C
// implementation in the same process.
//
// WHAT'S USED FROM LONGDUST (linked via FetchContent)
// ----------------------------------------------------
// From longdust.h (public API):
//   - ld_opt_t:       options struct (k-mer size, thresholds, etc.)
//   - ld_opt_init():  initialize options to defaults (k=7, T=0.6, ws=5000)
//   - ld_data_t:      opaque handle to internal scoring state
//   - ld_data_init(): allocate and precompute internal tables (incl. f[])
//   - ld_data_destroy(): free internal state
//
// From longdust.c (internal, accessed by struct replication):
//   - f[] array:      precomputed f(ℓ) table (Poisson expected value)
//   - opt->kmer:      k-mer size
//   - opt->ws:        window size (controls f[] table length)
//
// WHAT'S REIMPLEMENTED HERE
// -------------------------
// The following are reimplemented because longdust.c does not export them
// as standalone functions (they're embedded in the scanning algorithm):
//   - ld_test_nt4_table:      base-to-2bit encoding (copy of seq_nt4_table)
//   - ld_test_complement():   complement function (copy of seq_comp_tab)
//   - ld_test_q_one_strand(): single-strand Q(x)/ℓ computation
//   - ld_test_q_both_strands(): both-strand max (matching ld_dust2)
//
// IMPORTANT: The ld_test_data_s struct is a layout-compatible copy
// of longdust.c's internal ld_data_s struct. If longdust is updated,
// this struct MUST be kept in sync or the f[] pointer will be wrong.
// ============================================================================
#ifndef TESTS_BASE_LONGDUST_TEST_HELPERS_H_
#define TESTS_BASE_LONGDUST_TEST_HELPERS_H_

#ifdef __cplusplus
extern "C" {
#endif

// --- From longdust (linked into TestLancet2 via FetchContent) ---
#include "kdq.h"       // kdq_t — longdust's internal deque macro
#include "longdust.h"  // ld_opt_t, ld_data_t, ld_intv_t, public API

#include <math.h>
#include <stdlib.h>
#include <string.h>

// ============================================================================
// Replicated internal struct: ld_data_s (longdust.c lines 131-149)
//
// Longdust's ld_data_t is an opaque pointer (typedef struct ld_data_s
// ld_data_t) whose layout is defined only inside longdust.c. We replicate
// the struct layout here so we can access the precomputed f[] table.
//
// Field descriptions:
//   km         – kalloc memory pool (longdust's allocator)
//   opt        – pointer to the options (k-mer size, thresholds, etc.)
//   f          – precomputed f(ℓ) table; f[ℓ] = 4^k × f_single(ℓ/4^k)
//                This is the expected Σ log(c!) under the Poisson null.
//                Computed by ld_data_init() for ℓ = 0..ws+1.
//   c          – precomputed log(i) table; c[i] = log(i)
//                Used for incremental log-factorial: s += c[++ht[x]]
//   q          – deque for backward sweep positions (scanning algorithm)
//   ht         – hash table (flat array) for k-mer counts
//   max_test   – maximum count threshold for testing
//   n_for_pos  – number of forward positions stored
//   for_pos    – forward position/max array (scanning algorithm)
//   n_intv, m_intv – number/capacity of detected intervals
//   intv       – output interval array
// ============================================================================
KDQ_INIT(uint32_t)

typedef struct {
  int32_t pos;  // position in the scanned sequence
  double max;   // max Q/ℓ at this position
} ld_test_forpos_t;

struct ld_test_data_s {
  void* km;                   // kalloc memory pool
  ld_opt_t const* opt;        // pointer to scoring options
  double* f;                  // f[ℓ] = expected Σ log(c!) for ℓ k-mers (Poisson null)
  double* c;                  // c[i] = log(i), precomputed for i=0..max_count
  kdq_t(uint32_t) * q;        // deque used by backward/forward sweep
  uint16_t* ht;               // flat k-mer count array, size 4^k
  int32_t max_test;           // max count to trigger testing
  int32_t n_for_pos;          // number of forward pass entries
  ld_test_forpos_t* for_pos;  // forward pass position/score array
  int64_t n_intv;             // number of intervals found
  int64_t m_intv;             // allocated capacity for intervals
  ld_intv_t* intv;            // output array of low-complexity intervals
};

// ============================================================================
// Base encoding table — copied from longdust.c:97-114 (seq_nt4_table)
//
// Maps ASCII characters to 2-bit DNA codes:
//   'A'/'a' → 0    (binary 00)
//   'C'/'c' → 1    (binary 01)
//   'G'/'g' → 2    (binary 10)
//   'T'/'t' → 3    (binary 11)
//   everything else → 4  (sentinel: not a valid DNA base)
//
// The sentinel value 4 is used to detect N bases and other invalid
// characters, which reset the k-mer sliding window.
// ============================================================================
static unsigned char ld_test_nt4_table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x00..0x0F
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x10..0x1F
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x20..0x2F (space, punctuation)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x30..0x3F (digits)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x40..0x4F (@,A,B,C,...,O) — A=0,C=1,G=2
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x50..0x5F (P,...,T,...,_)   — T=3
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x60..0x6F (`,a,b,c,...,o) — a=0,c=1,g=2
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x70..0x7F (p,...,t,...,DEL) — t=3
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0x80..0xFF: all invalid
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

// ============================================================================
// Complement function — reimplements longdust.c:116-125 (seq_comp_tab)
//
// Returns the Watson-Crick complement of a DNA base, preserving case.
// Non-DNA characters (including N) return 'N'.
// ============================================================================
inline static char ld_test_complement(char c) {
  switch (c) {
    case 'A':
      return 'T';
    case 'a':
      return 't';
    case 'T':
      return 'A';
    case 't':
      return 'a';
    case 'C':
      return 'G';
    case 'c':
      return 'g';
    case 'G':
      return 'C';
    case 'g':
      return 'c';
    default:
      return 'N';
  }
}

// ============================================================================
// ld_test_q_one_strand: compute Q(x)/ℓ on a single strand
//
// This reimplements the core scoring logic that longdust uses internally
// during its backward/forward sweep (longdust.c:220-240). Longdust does
// this incrementally as it slides a window; we do it in batch on a fixed
// sequence. The math is identical:
//
//   Q(x) = Σ_t log(c_x(t)!) − f(ℓ)
//   q(x) = max(0, Q(x) / ℓ)
//
// Variables (using longdust's naming conventions):
//   d       – cast to our replicated struct to access internal fields
//   k       – k-mer length (d->opt->kmer, typically 7)
//   mask    – bitmask to keep only the low 2k bits: (1 << 2k) - 1
//   n_kmer  – total possible k-mers: 4^k = 1 << 2k (16,384 for k=7)
//   ht      – flat count array indexed by k-mer code (longdust calls it ht)
//   x       – current k-mer as a 2k-bit integer, built by shift+OR
//   l       – run length: counts consecutive valid bases without N
//   valid   – total valid k-mers extracted (this is ℓ in the formula)
//   b       – encoded base (0-3 for ACGT, 4 for invalid/N)
//   sum_log_fact – accumulates Σ log(c!)
//   f_val   – f(ℓ) from longdust's precomputed table
//   q       – raw score Q(x) = sum_log_fact - f_val
// ============================================================================
inline static double ld_test_q_one_strand(ld_data_t const* ld, char const* seq, int seq_len) {
  // Cast opaque ld_data_t to our replicated struct to access f[] table
  const struct ld_test_data_s* d = (const struct ld_test_data_s*)ld;

  // Read k-mer parameters from longdust's options
  int const k = d->opt->kmer;                 // k-mer size (typically 7)
  uint32_t const mask = (1U << (2 * k)) - 1;  // bitmask for 2k-bit k-mer
  uint32_t const n_kmer = 1U << (2 * k);      // number of possible k-mers (4^k)

  // Allocate k-mer count array (zeroed), same as longdust's ht[] array
  uint16_t* ht = (uint16_t*)calloc(n_kmer, sizeof(uint16_t));

  uint32_t x = 0;  // current k-mer as a 2k-bit integer
  int l = 0;       // run length of consecutive valid bases
  int valid = 0;   // total number of valid k-mers extracted (= ℓ)

  // ── Step 1: Slide through sequence and count k-mers ──────────────────
  // This matches longdust.c:215-228 (the forward scan in ld_dust1)
  for (int i = 0; i < seq_len; i++) {
    int b = ld_test_nt4_table[(unsigned char)seq[i]];  // encode base
    if (b < 4) {
      // Valid base: shift into k-mer integer and mask to k positions
      x = (x << 2 | b) & mask;
      // Once we've seen k consecutive valid bases, we have a full k-mer
      if (++l >= k) {
        ht[x]++;  // increment count for this k-mer
        valid++;  // increment total valid k-mer count
      }
    } else {
      l = 0;  // N or invalid base: reset run, next k-mer needs k new bases
    }
  }

  // No valid k-mers → score is 0 (sequence too short or all N's)
  if (valid == 0) {
    free(ht);
    return 0.0;
  }

  // ── Step 2: Compute Σ_t log(c_x(t)!) ────────────────────────────────
  // Longdust accumulates this incrementally as: s += ld->c[++ht[x>>1]]
  // where c[i] = log(i). That's equivalent to summing log(2)+log(3)+...+log(c)
  // = log(c!) for each k-mer with count c ≥ 2.
  //
  // We do the same thing in batch: for each k-mer, if count ≥ 2,
  // add log(2) + log(3) + ... + log(count) = log(count!).
  double sum_log_fact = 0.0;
  for (uint32_t t = 0; t < n_kmer; t++) {
    if (ht[t] >= 2) {
      for (int j = 2; j <= ht[t]; j++) {
        sum_log_fact += log((double)j);  // Σ log(j) = log(count!)
      }
    }
  }

  // ── Step 3: Subtract f(ℓ) and normalize ──────────────────────────────
  // f(ℓ) is the expected value of Σ log(c!) under the random (Poisson) model.
  // Longdust precomputes f[] for ℓ = 0..ws+1 in ld_data_init().
  // The raw score is: Q = sum_log_fact - f(ℓ)
  // The per-k-mer score is: q = max(0, Q / ℓ)
  double f_val = (valid <= d->opt->ws + 1) ? d->f[valid] : 0.0;
  double q = sum_log_fact - f_val;

  free(ht);
  return (valid > 0 && q > 0.0) ? q / valid : 0.0;
}

// ============================================================================
// ld_test_q_both_strands: compute max(q_fwd, q_rev) matching ld_dust2
//
// This reimplements the both-strand logic from longdust.c:361-411 (ld_dust2).
// Longdust runs ld_dust1 on the forward strand, then builds the reverse
// complement and runs ld_dust1 again, merging the resulting intervals.
// We just compute the raw score on each strand and take the maximum.
// ============================================================================
inline static double ld_test_q_both_strands(ld_data_t const* ld, char const* seq, int seq_len) {
  // Score the forward strand
  double score_fwd = ld_test_q_one_strand(ld, seq, seq_len);

  // Build reverse complement — matches longdust.c:374-376
  // longdust reverses the sequence and complements each base
  char* rev = (char*)malloc(seq_len + 1);
  for (int i = 0; i < seq_len; i++) {
    rev[seq_len - 1 - i] = ld_test_complement(seq[i]);
  }
  rev[seq_len] = '\0';

  // Score the reverse complement strand
  double score_rev = ld_test_q_one_strand(ld, rev, seq_len);
  free(rev);

  // Return the higher of the two scores (matching ld_dust2 behavior)
  return score_fwd > score_rev ? score_fwd : score_rev;
}

#ifdef __cplusplus
}
#endif

#endif  // TESTS_BASE_LONGDUST_TEST_HELPERS_H_
