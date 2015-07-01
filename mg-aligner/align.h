//-------------------------------------------------------
// Common Alignment Functionality
// Victoria Popic (viq@stanford.edu), 2 Apr 2012
//-------------------------------------------------------

#ifndef BWTSNP_ALIGN_H
#define BWTSNP_ALIGN_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

#define READ_BATCH_SIZE 0x40000

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define MAX_DIFF 5
#define ALN_PATH_ALLOC 256
#define MAX_SNPS 5

#define ALN_NOMATCH 0
#define ALN_UNIQUE  1
#define ALN_REPEAT  2

#define MAPQ_CONFIDENT 10

#define NUM_PRECALC 16777216
#define PRECALC_INTERVAL_LENGTH 12

// SA interval
typedef struct sa_intv_ {
	bwtint_t L; // lower bound
	bwtint_t U; // upper bound
	struct sa_intv_* next_intv;
} sa_intv_t;

// list of SA intervals
// invariant: intervals are disjoint and sorted
typedef struct {
	int size;
	sa_intv_t* first_intv;
	sa_intv_t* last_intv;
} sa_intv_list_t;

typedef struct {
	int max_diff;
	int max_gapo;
	int max_gape;
	int max_entries;
	// scoring parameters
	int mm_score; // mismatch
	int gapo_score; // open gap
	int gape_score; // extend gap

	// length of the seed (seed := first seed_length chars of the read)
	int seed_length;
	// maximum number of differences in the seed
	int max_diff_seed;
	// minimum number of required top hits to not proceed with suboptimal alignments
	int max_best;
	// bp length at each end where indels will not be allowed
	int no_indel_length;

	// number of contiguous N's allowed to match
	int matched_Ncontig;

	// flag indicating whether to use pre-computed intervals
	int use_precalc;

	// flag indicating alignment with multi-reference
	int is_multiref;

	// multi-threading
	int n_threads;

} aln_params_t;

typedef struct {
	int score;
	bwtint_t L, U;
	int num_mm;
	int num_gapo;
	int num_gape;
	int num_snps;
	int aln_length;
	char* aln_path;
} aln_t;

typedef struct {
	// number of alignment entries
	int num_entries;
	// number of allocated alignment entries
	int max_entries;
	aln_t* entries;
} alns_t;

typedef struct {
	bwtint_t L;
	bwtint_t U; // (L,U): SA interval of [i,n-1]
	uint32_t num_mm:8, num_gapo:8, num_gape:8, num_snps:8;
	uint32_t score:8, i:8, state:8, aln_length:8;

	//uint8_t score; // aln score so far
	//uint8_t i; // position in the read
	//uint8_t num_mm;
	//uint8_t num_gapo;
	//uint8_t num_gape;
	//uint8_t state;
	//int n_seed_mm;
	//uint8_t last_diff_pos;
	//uint8_t padding;

	// edit transcript
	//uint8_t aln_length;
	char aln_path[ALN_PATH_ALLOC];
} aln_entry_t;


// Read Alignment
int align_reads(char* fastaFname, char* readsFname, char* alnsFname, aln_params_t* params);

// SA Intervals Operations
void add_sa_interval(sa_intv_list_t* intv_list, bwtint_t L, bwtint_t U);
void clear_sa_interval_list(sa_intv_list_t* intv_list);
void free_sa_interval_list(sa_intv_list_t* intv_list);
void print_sa_interval_list(sa_intv_list_t* intv_list);
void store_sa_interval_list(sa_intv_list_t* intv_list, FILE* saFile);
void load_sa_interval_list(sa_intv_list_t* intv_list, FILE* saFile);
int read2index(char* read, int readLen);

// Alignments Operations
void set_default_aln_params(aln_params_t* params);
alns_t* init_alignments();
alns_t* sa_intervals2alns(sa_intv_list_t* intv_list, int aln_length);
void free_alignments(alns_t* alns);
void reset_alignments(alns_t* alns);
void add_alignment(aln_entry_t* e, bwtint_t L, bwtint_t U, int score, alns_t* alns, const aln_params_t* params);
void print_alignments(alns_t* alns);
void alns2alnf(alns_t* alns, FILE* alnFile);
void alns2alnf_bin(alns_t* alns, FILE* alnFile);
alns_t* alnsf2alns(int* num_alns, char *alnFname);
alns_t* alnsf2alns_bin(int* n_alns, char *alnFname);
void eval_alns(char *fastaFname, char *readsFname, char *alnFname, int is_multiref, int max_diff);
void alns2sam(char *fastaFname, char *readsFname, char *alnsFname, char* samFname, int is_multiref, int max_diff);

#endif
