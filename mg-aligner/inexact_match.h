//-------------------------------------------------------
// Inexact Matching
// Victoria Popic (viq@stanford.edu), Apr 2012
//-------------------------------------------------------

#ifndef BWTSNP_INEXACT_MATCH_H
#define BWTSNP_INEXACT_MATCH_H

#include <stdint.h>

typedef struct {
	int num_diff;
	int sa_intv_width;
} diff_lower_bound_t;

// stores alignment entries with the same score
typedef struct {
	// number of alignment entries in the bucket
	int num_entries;
	// number of allocated entries
	int max_entries;
	aln_entry_t* entries;
} heap_bucket_t;

typedef struct {
	// index into the best score non-empty bucket
	int best_score;
	// number of different possible alignment scores
	int num_buckets;
	// total number of alignment entries in all the buckets
	int num_entries;
	// store alignment entries with the same score
	heap_bucket_t* buckets;
} priority_heap_t;

/* Operations */

// Inexact matching operations
int align_reads_inexact(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals, aln_params_t* params, char* alnsFname);
int align_reads_inexact_parallel(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals_table, aln_params_t* params, char* alnFname);
void inexact_match(bwt_t * BWT, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, const aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns);

// Priority heap operations
priority_heap_t *heap_init(const aln_params_t *p);
void heap_free(priority_heap_t *heap);
void heap_reset(priority_heap_t *heap);
void heap_push(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U, const int num_mm, const int num_gapo, const int num_gape,
		const int state, const int is_diff, const int aln_length, const char* aln_path, const aln_params_t *params);
void heap_pop(priority_heap_t *heap, aln_entry_t* e);

#endif
