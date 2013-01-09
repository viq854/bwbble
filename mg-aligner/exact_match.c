//-------------------------------------------------------
// Exact Matching
// Victoria Popic (viq@stanford.edu), Apr 2012
//-------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "align.h"
#include "bwt.h"
#include "exact_match.h"

int align_reads_exact(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals, const aln_params_t* params, char* alnFname) {
	printf("BWT-SNP Exact Alignment...\n");
	sa_intv_list_t** sa_intervals = (sa_intv_list_t**) malloc(reads->count*sizeof(sa_intv_list_t*));
	clock_t t = clock();
	for(int i = 0; i < reads->count; i++) {
		// match with the forward and reverse complement references
		if(params->use_precalc) {
			exact_match_precalc(BWT, &reads->reads[i], precalc_sa_intervals, &sa_intervals[i], params);
		} else {
			exact_match(BWT, &reads->reads[i], &sa_intervals[i], params);
		}
	}
	printf("Exact matching time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	t = clock();
	alns_t** alns = (alns_t**) malloc(reads->count*sizeof(alns_t*));
	for(int i = 0; i < reads->count; i++) {
		reads->reads[i].alns = sa_intervals2alns(sa_intervals[i], reads->reads[i].len);
		free_sa_interval_list(sa_intervals[i]);
	}
	// write the results to file
	FILE* alnFile = (FILE*) fopen(alnFname, "a+");
	if (alnFile == NULL) {
		printf("align_reads_exact: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}
	for(int i = 0; i < reads->count; i++) {
		read_t* read = &reads->reads[i];
		alns2alnf(read->alns, alnFile);
		free_alignments(read->alns);
	}
	fclose(alnFile);
	printf("Storing results time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	free(sa_intervals);
	free(alns);
	return 0;
}

// returns 1 if any matches were found and the list of matching SA intervals
int exact_match(bwt_t *BWT, read_t* read, sa_intv_list_t** sa_intervals, const aln_params_t* params) {
	return exact_match_bounded(BWT, read->seq, read->len, 0, BWT->length-1, read->len-1, sa_intervals, params);
}

int exact_match_1to1_bounded(bwt_t *BWT, const char* read, const int readLen, const bwtint_t l, const bwtint_t u, const int i, bwtint_t *sa_begin, bwtint_t *sa_end);

// returns 1 if any matches were found with read[0,i] starting with the SA interval (l,u)
// returns the list of matching SA intervals
int exact_match_bounded(bwt_t *BWT, char* read, int readLen, bwtint_t l, bwtint_t u, int i, sa_intv_list_t** sa_intervals, const aln_params_t* params) {
	sa_intv_list_t* intv_list_curr = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));

	if(!params->is_multiref) {
		bwtint_t L, U;
		int matched = exact_match_1to1_bounded(BWT, read, readLen, l, u, i, &L, &U);
		if(matched > 0) {
			add_sa_interval(intv_list_curr, L, U);
		}
		*sa_intervals = intv_list_curr;
		return (intv_list_curr->size != 0);
	}

	sa_intv_list_t* intv_list_next = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
	add_sa_interval(intv_list_curr, l, u);

	for (int r = i; r >= 0; r--) {
		unsigned char c = read[r];
		if(c == nt4_table[(int) 'N']) {
			clear_sa_interval_list(intv_list_curr);
			break; // N in the read is considered a mismatch
		}
		sa_intv_t* intv = intv_list_curr->first_intv;
		for(int s = 0; s < intv_list_curr->size; s++) {
			if(params->is_multiref) {
				for(int b = 0; b < BASES_PER_NUCLEOTIDE; b++) {
					unsigned char base = nucl_bases_table[c][b];
					if(base == 10) continue; // do not match with N's in the reference
					bwtint_t L = BWT->C[base] + O(BWT, base, intv->L-1) + 1;
					bwtint_t U = BWT->C[base] + O(BWT, base, intv->U);
					if (L <= U) {
						add_sa_interval(intv_list_next, L, U);
					}
				}
			} else {
				unsigned char base = nt4_gray[c];
				bwtint_t L = BWT->C[base] + O(BWT, base, intv->L-1) + 1;
				bwtint_t U = BWT->C[base] + O(BWT, base, intv->U);
				if (L <= U) {
					add_sa_interval(intv_list_next, L, U);
				}
			}
			intv = intv->next_intv;
		}
		sa_intv_list_t* tmp = intv_list_curr;
		intv_list_curr = intv_list_next;
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);
		if(intv_list_curr->size == 0) break;
	}
	*sa_intervals = intv_list_curr;
	free(intv_list_next);
	return (intv_list_curr->size != 0);
}

// returns 1 if any matches were found and the list of matching SA intervals
int exact_match_precalc(bwt_t *BWT, read_t* read, sa_intv_list_t* intv_table, sa_intv_list_t** sa_intervals, const aln_params_t* params) {
	if(read->len < PRECALC_INTERVAL_LENGTH) {
		return exact_match_bounded(BWT, read->seq, read->len, 0, BWT->length-1, read->len-1, sa_intervals, params);
	}

	int read_index = read2index(read->seq, read->len);
	sa_intv_list_t* intv_list_precalc = &(intv_table[read_index]);
	sa_intv_list_t* intv_list_curr = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
	sa_intv_list_t* intv_list_next = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));

	if (intv_list_precalc->size != 0) {
		sa_intv_t* intv = intv_list_precalc->first_intv;
		for(int i = 0; i < intv_list_precalc->size; i++) {
			add_sa_interval(intv_list_curr, intv->L, intv->U);
			intv = intv->next_intv;
		}
	}

	for (int r = read->len - PRECALC_INTERVAL_LENGTH - 1; r >= 0; r--) {
		unsigned char c = read->seq[r];
		if(c == nt4_table[(int) 'N']) {
			clear_sa_interval_list(intv_list_curr);
			break; // N in the read is considered a mismatch
		}
		sa_intv_t* intv = intv_list_curr->first_intv;
		for(int s = 0; s < intv_list_curr->size; s++) {
			if(params->is_multiref) {
				for(int b = 0; b < BASES_PER_NUCLEOTIDE; b++) {
					unsigned char base = nucl_bases_table[c][b];
					if(base == 10) continue; // do not match with N's in the reference
					bwtint_t L = BWT->C[base] + O(BWT, base, intv->L-1) + 1;
					bwtint_t U = BWT->C[base] + O(BWT, base, intv->U);
					if (L <= U) {
						add_sa_interval(intv_list_next, L, U);
					}
				}
			} else {
				unsigned char base = nt4_gray[c];
				bwtint_t L = BWT->C[base] + O(BWT, base, intv->L-1) + 1;
				bwtint_t U = BWT->C[base] + O(BWT, base, intv->U);
				if (L <= U) {
					add_sa_interval(intv_list_next, L, U);
				}
			}
			intv = intv->next_intv;
		}
		sa_intv_list_t* tmp = intv_list_curr;
		intv_list_curr = intv_list_next;
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);
		if(intv_list_curr->size == 0) break;
	}
	*sa_intervals = intv_list_curr;
	free(intv_list_next);
	return (intv_list_curr->size != 0);
}

// 1-to-1 matching
// returns the number of matches and the sa interval
int exact_match_1to1(bwt_t *BWT, read_t* read, bwtint_t *sa_begin, bwtint_t *sa_end) {
	bwtint_t L = 0;
	bwtint_t U = BWT->length;
	for (int i = (read->len - 1); i >= 0; i--) {
		unsigned char c = read->seq[i];
		L = BWT->C[c] + O(BWT, c, L-1) + 1;
		U = BWT->C[c] + O(BWT, c, U);
		if (L > U) break;
	}
	if (L > U) return 0; // no match
	(*sa_begin) = L;
	(*sa_end) = U;
	return U - L + 1;
}

int exact_match_1to1_bounded(bwt_t *BWT, const char* read, const int readLen, const bwtint_t l, const bwtint_t u, const int i, bwtint_t *sa_begin, bwtint_t *sa_end) {
	bwtint_t L = l;
	bwtint_t U = u;
	for (int j = i; j >= 0; j--) {
		if(read[j] > 3) { // N in the read is considered a mismatch
			return 0;
		}
		const unsigned char c = nt4_gray[(int) read[j]];
		bwtint_t occL, occU;
		if((L-1) == U) {
			occL = O(BWT, c, L-1);
			occU = occL;
		} else {
			occL = O(BWT, c, L-1);
			occU = O(BWT, c, U);
			//O_LU(BWT, c, L-1, U, &occL, &occU);
		}
		L = BWT->C[c] + occL + 1;
		U = BWT->C[c] + occU;
		if (L > U) {
			return 0; // no match
		}
	}
	(*sa_begin) = L;
	(*sa_end) = U;
	return U - L + 1;
}
