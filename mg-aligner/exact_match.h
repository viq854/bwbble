//-------------------------------------------------------
// Exact Matching
// Victoria Popic (viq@stanford.edu), Apr 2012
//-------------------------------------------------------

#ifndef BWTSNP_EXACT_MATCH_H
#define BWTSNP_EXACT_MATCH_H

#include <stdint.h>

// Exact matching (using backwards search)
int align_reads_exact(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals, const aln_params_t* params, char* alnsFname);
int exact_match(bwt_t *BWT, read_t* read, sa_intv_list_t** sa_intervals, const aln_params_t* params);
int exact_match_bounded(bwt_t *BWT, char* read, int readLen, bwtint_t L, bwtint_t U, int i, sa_intv_list_t** sa_intervals, const aln_params_t* params);
int exact_match_precalc(bwt_t *BWT, read_t* read, sa_intv_list_t* intv_table, sa_intv_list_t** sa_intervals, const aln_params_t* params);
int exact_match_1to1(bwt_t *BWT, read_t* read, bwtint_t *sa_begin, bwtint_t *sa_end);

#endif
