//-------------------------------------------------------
// Inexact Matching
// Victoria Popic (viq@stanford.edu), Apr 2012
//-------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "bwt.h"
#include "align.h"
#include "inexact_match.h"
#include "exact_match.h"
#include <omp.h>

void calculate_d(bwt_t* BWT, char* read, const int readLen, diff_lower_bound_t* D, aln_params_t* params);

static inline int aln_score(const int m, const int o, const int e, const aln_params_t *p) {
	return m*p->mm_score + o*p->gapo_score + e*p->gape_score;
}

int align_reads_inexact(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals_table, aln_params_t* params, char* alnFname) {
	printf("BWBBLE Inexact Alignment...\n");
	FILE* alnFile = (FILE*) fopen(alnFname, "a+b");
	if (alnFile == NULL) {
		printf("align_reads_inexact: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}
	// lower bound on the number of differences at each position in the read
	diff_lower_bound_t* D = (diff_lower_bound_t*) calloc(reads->max_len+1, sizeof(diff_lower_bound_t));
	// lower bound for the read seed positions
	diff_lower_bound_t* D_seed = (diff_lower_bound_t*) calloc(params->seed_length+1, sizeof(diff_lower_bound_t));

	priority_heap_t* heap = heap_init(params);

	// process the reads in batches
	int num_processed = 0;
	while(num_processed < reads->count) {
		clock_t t = clock();
		int batch_size = ((reads->count - num_processed) > READ_BATCH_SIZE ) ? READ_BATCH_SIZE : (reads->count - num_processed);

		for(int i = num_processed; i < num_processed + batch_size; i++) {
			read_t* read = &reads->reads[i];
			read->alns = init_alignments();
			sa_intv_list_t* precalc_sa_intervals = NULL;
			if(params->use_precalc) {
				// discard reads that have N's in the last PRECALC_INTERVAL_LENGTH bases (<0 result from read_index)
				int read_index = read2index(read->rc, read->len);
				if(read_index < 0) {
					continue;
				}
				precalc_sa_intervals = &(precalc_sa_intervals_table[read_index]);
			}

			// align read with forward reference <=> read reverse complement with BWT reverse complement
			// align read with reverse complement reference <=> read reverse complement with BWT forward
			calculate_d(BWT, read->seq, read->len, D, params);
			if(params->seed_length && read->len > params->seed_length) {
				calculate_d(BWT, read->seq, params->seed_length, D_seed, params);
			}
			inexact_match(BWT, read->rc, read->len, heap, precalc_sa_intervals, params, D, D_seed, read->alns);
		}
		printf("Processed %d reads. Inexact matching time: %.2f sec.", num_processed+batch_size, (float)(clock() - t) / CLOCKS_PER_SEC);

		// write the results to file
		clock_t ts = clock();
		for(int i = num_processed; i < num_processed + batch_size; i++) {
			read_t* read = &reads->reads[i];
			alns2alnf_bin(read->alns, alnFile);
			free_alignments(read->alns);
			free(read->seq);
			free(read->rc);
			free(read->qual);
			read->seq = read->rc = read->qual = NULL;
		}
		printf("Storing results time: %.2f sec\n", (float)(clock() - ts) / CLOCKS_PER_SEC);
		num_processed += batch_size;
	}

	free(D);
	free(D_seed);
	heap_free(heap);
	fclose(alnFile);
	return 0;
}

// naive parallelization scheme (1 thread <=> 1 read)
int align_reads_inexact_parallel(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals_table, aln_params_t* params, char* alnFname) {
	printf("BWT-SNP Inexact Alignment...\n");
	FILE* alnFile = (FILE*) fopen(alnFname, "a+");
	if (alnFile == NULL) {
		printf("align_reads_inexact: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}

	// process the reads in batches
	int num_processed = 0;
	while(num_processed < reads->count) {
		clock_t t = clock();
		int batch_size = ((reads->count - num_processed) > READ_BATCH_SIZE ) ? READ_BATCH_SIZE : (reads->count - num_processed);

		omp_set_num_threads(params->n_threads);
		int tid, n_threads, chunk_start, chunk_end;
		diff_lower_bound_t* D, * D_seed;
		priority_heap_t* heap;
		#pragma omp parallel private(tid, n_threads, chunk_start, chunk_end, D, D_seed, heap)
		{
			tid = omp_get_thread_num();
			n_threads = omp_get_num_threads();
			chunk_start = tid * batch_size / n_threads;
			chunk_end = (tid + 1) * batch_size / n_threads;

			// lower bound on the number of differences at each position in the read
			D = (diff_lower_bound_t*) calloc(reads->max_len+1, sizeof(diff_lower_bound_t));
			// lower bound for the read seed positions
			D_seed = (diff_lower_bound_t*) calloc(params->seed_length+1, sizeof(diff_lower_bound_t));
			// partial alignments min-heap
			heap = heap_init(params);

			for (int i = num_processed + chunk_start; i < num_processed + chunk_end; i++) {
				read_t* read = &reads->reads[i];
				read->alns = init_alignments();
				sa_intv_list_t* precalc_sa_intervals = NULL;
				if(params->use_precalc) {
					// discard reads that have N's in the last PRECALC_INTERVAL_LENGTH bases (<0 result from read_index)
					int read_index = read2index(read->rc, read->len);
					if(read_index < 0) {
						continue;
					}
					precalc_sa_intervals = &(precalc_sa_intervals_table[read_index]);
				}

				// align read with forward reference <=> read reverse complement with BWT reverse complement
				// align read with reverse complement reference <=> read reverse complement with BWT forward
				calculate_d(BWT, read->seq, read->len, D, params);
				if(params->seed_length && read->len > params->seed_length) {
					calculate_d(BWT, read->seq, params->seed_length, D_seed, params);
				}
				inexact_match(BWT, read->rc, read->len, heap, precalc_sa_intervals, params, D, D_seed, read->alns);
			}
			free(D);
			free(D_seed);
			heap_free(heap);
		}
		printf("Processed %d reads. Inexact matching time: %.2f sec.", num_processed+batch_size, (float)(clock() - t) / CLOCKS_PER_SEC);

		// write the results to file
		clock_t ts = clock();
		for(int i = num_processed; i < num_processed + batch_size; i++) {
			read_t* read = &reads->reads[i];
			alns2alnf_bin(read->alns, alnFile);
			free_alignments(read->alns);
			free(read->seq);
			free(read->rc);
			free(read->qual);
			read->seq = read->rc = read->qual = NULL;
		}
		printf("Storing results time: %.2f sec\n", (float)(clock() - ts) / CLOCKS_PER_SEC);
		num_processed += batch_size;
	}
	fclose(alnFile);
	return 0;
}

// Lower bound on the number of differences at each position in the read (used in BWA)
void calculate_d(bwt_t *BWT, char* read, int readLen, diff_lower_bound_t* D, aln_params_t* params) {
	int z = 0;
	bwtint_t L = 0;
	bwtint_t U = BWT->length-1;

	if(!params->is_multiref) {
		for (int i = readLen-1; i >= 0; i--) {
			unsigned char c = nt4_gray[(int) read[i]];
			if(c == 10) {
				L = 0;
				U = BWT->length-1;
				z++;
			} else {
				bwtint_t occL, occU;
				if((L-1) == U) {
					occL = O(BWT, c, L-1);
					occU = occL;
				} else {
					occL = O(BWT, c, L-1);
					occU = O(BWT, c, U);
				}
				L = BWT->C[c] + occL + 1;
				U = BWT->C[c] + occU;
				if(L > U) {
					L = 0;
					U = BWT->length-1;
					z++;
				}
			}
			D[readLen-1-i].num_diff = z;
			D[readLen-1-i].sa_intv_width = U-L+1;
		}
		D[readLen].sa_intv_width = 0;
		D[readLen].num_diff = ++z;
		return;
	}

	sa_intv_list_t* intv_list_curr = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
	sa_intv_list_t* intv_list_next = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
	add_sa_interval(intv_list_curr, L, U);

	// backwards search
	for (int i = readLen-1; i >= 0; i--) {
		unsigned char c = read[i];
		int num_matches = 0;
		if(c > 3) {
			clear_sa_interval_list(intv_list_curr); // N in the read is considered a mismatch
		} else {
			sa_intv_t* intv = intv_list_curr->first_intv;
			for(int s = 0; s < intv_list_curr->size; s++) {
				for(int b = 0; b < BASES_PER_NUCLEOTIDE; b++) {
					unsigned char base = nucl_bases_table[c][b];
					if(base == 10) continue; // do not match with N's in the reference
					bwtint_t L = BWT->C[base] + O(BWT, base, intv->L-1) + 1;
					bwtint_t U = BWT->C[base] + O(BWT, base, intv->U);
					if (L <= U) {
						num_matches += U - L + 1;
						add_sa_interval(intv_list_next, L, U);
					}
				}
				intv = intv->next_intv;
			}
		}
		sa_intv_list_t* tmp = intv_list_curr;
		intv_list_curr = intv_list_next;
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);

		// no matches
		if(intv_list_curr->size == 0) {
			add_sa_interval(intv_list_curr, 0, BWT->length-1);
			z++;
			num_matches = U - L + 1;
		}
		D[readLen-1-i].num_diff = z;
		D[readLen-1-i].sa_intv_width = num_matches;
	}

	D[readLen].sa_intv_width = 0;
	D[readLen].num_diff = ++z;

	free_sa_interval_list(intv_list_curr);
	free_sa_interval_list(intv_list_next);
}

void inexact_match(bwt_t * BWT, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals,
		const aln_params_t *params, diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns) {

	// discard any reads where the number of N's is higher than the number of allowed mismatches
	int countN = 0;
	for (int i = 0; i < readLen; i++) {
		if (read[i] > 3) countN++;
	}
	if (countN > params->max_diff) {
		return;
	}

	heap_reset(heap);
	if(precalc_sa_intervals != NULL) {
		if (precalc_sa_intervals->size != 0) {
			char aln_path[ALN_PATH_ALLOC] = { 0 };
            sa_intv_t* intv = precalc_sa_intervals->first_intv;
			for(int i = 0; i < precalc_sa_intervals->size; i++) {
				heap_push(heap, readLen - PRECALC_INTERVAL_LENGTH, intv->L, intv->U, 0, 0, 0, 0, 0, PRECALC_INTERVAL_LENGTH-1, aln_path, params);
				intv = intv->next_intv;
			}
		} else {
        	return; // no matches
       	}
	} else {
		heap_push(heap, readLen, 0, BWT->length-1, 0, 0, 0, 0, 0, 0, 0, params);
	}
	
	int best_score = aln_score(params->max_diff+1, params->max_gapo+1, params->max_gape+1, params); 
	int best_diff = params->max_diff + 1;
	int max_diff = params->max_diff;
	int num_best = 0;
	int max_entries = 0;

	int total_entries = 0;
	int last_num_entries = 1;

	while (heap->num_entries != 0) {

		total_entries += heap->num_entries - last_num_entries + 1;
		last_num_entries = heap->num_entries;

		if(heap->num_entries > max_entries) max_entries = heap->num_entries;
		if (heap->num_entries > params->max_entries) {
			break;
		}

		// pop best entry
		aln_entry_t e_;
		heap_pop(heap, &e_);
		aln_entry_t* e = &e_;

		// case 1 optimization
		if(e->score > (best_score + params->mm_score)) {
			break;
		}
		int diff_left = max_diff - e->num_mm - e->num_gapo - e->num_gape;
		if (diff_left < 0) {
			continue;
		}
		// apply the lower bound estimate
		if ((e->i > 0) && (diff_left < D[e->i-1].num_diff)) {
			continue;
		}
		// apply the seed
		int diff_left_seed = params->max_diff_seed - e->num_mm - e->num_gapo - e->num_gape;
		int seed_index = e->i - (readLen - params->seed_length);
		//if ((seed_index == 0) && (diff_left_seed < 0)) {
			//continue;
		//}
		if ((seed_index > 0) && (diff_left_seed < D_seed[seed_index-1].num_diff)) {
			continue;
		}

		// check if this entry is a hit
		if (e->i == 0) { // all characters have been matched
			int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
			if(alns->num_entries == 0) { // first hit has the best score
				best_score = score;
				best_diff = e->num_mm + e->num_gapo + e->num_gape;
				max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
			}
			if(score == best_score) {
				num_best += e->U - e->L + 1;
			} else if(num_best > params->max_best) {
				break;
			}
			add_alignment(e, e->L, e->U, score, alns, params);
			continue;
		} else if (diff_left == 0) { // check if the remaining characters are an exact match
			sa_intv_list_t* sa_intervals;
			if (exact_match_bounded(BWT, read, readLen, e->L, e->U, e->i-1, &sa_intervals, params) > 0) {
				int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
				if(alns->num_entries == 0) { // first hit has the best score
					best_score = score;
					best_diff = e->num_mm + e->num_gapo + e->num_gape;
					max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
				}
				if(score == best_score) {
					sa_intv_t* intv = sa_intervals->first_intv;
					for(int k = 0; k < sa_intervals->size; k++) {
						num_best += intv->U - intv->L + 1;
						intv = intv->next_intv;
					}
				} else if(num_best > params->max_best) {
					break;
				}

				// add the matches to the path (don't need to change the path since STATE_M=0)
				e->aln_length += e->i;
				sa_intv_t* intv = sa_intervals->first_intv;
				for(int k = 0; k < sa_intervals->size; k++) {
					add_alignment(e, intv->L, intv->U, score, alns, params);
					intv = intv->next_intv;
				}
			}

			free_sa_interval_list(sa_intervals);
			continue;
		}

		bwtint_t L[ALPHABET_SIZE] = { 0 };
		bwtint_t U[ALPHABET_SIZE] = { 0 };
		int alphabet_size = ALPHABET_SIZE;
		int is_multiref = params->is_multiref;
		if(params->is_multiref) {
			O_alphabet(BWT, e->L-1, ALPHABET_SIZE, L, 1);
			O_alphabet(BWT, e->U, ALPHABET_SIZE, U, 0);
		} else {
			O_actg_alphabet(BWT, e->L-1, L, 1);
			O_actg_alphabet(BWT, e->U, U, 0);
			alphabet_size = NUM_NUCLEOTIDES+1; // (to be able to skip j=0 bellow)
			is_multiref = 0;
		}

		// check if any more differences can be introduced (includes BWA heuristics)
		int allow_diff = 1;
		int allow_indels = 1;
		int allow_mm = 1;
		int allow_open = 1;
		int allow_extend = 1;

		// check lower bound
		if(e->i-1 > 0) {
			if((diff_left - 1) < D[e->i-2].num_diff) {
				allow_diff = 0;
			} else if(((D[e->i-1].num_diff == diff_left - 1) && (D[e->i-2].num_diff == diff_left - 1))
					&& (D[e->i-1].sa_intv_width == D[e->i-2].sa_intv_width)) {
				allow_mm = 0;
			}
		}
		// check seed
		if(seed_index-1 > 0) {
			if((diff_left_seed - 1) < D_seed[seed_index-2].num_diff) {
				allow_diff = 0;
			} else if((D_seed[seed_index-1].num_diff == diff_left_seed - 1) && (D_seed[seed_index-2].num_diff == diff_left_seed -1)
					&& (D_seed[seed_index-1].sa_intv_width == D_seed[seed_index-2].sa_intv_width)) {
				allow_mm = 0;
			}
		}

		// check indels
		int tmp = e->num_gapo + e->num_gape;
		if((e->i-1 < (params->no_indel_length + tmp)) || ((readLen - (e->i-1)) < (params->no_indel_length + tmp))) {
			allow_indels = 0;
		}
		if((e->num_gapo >= params->max_gapo) && (e->num_gape >= params->max_gape)) {
			allow_indels = 0;
		}
		if((e->num_gapo >= params->max_gapo)) {
			allow_open = 0;
		}
		if((e->num_gape >= params->max_gape)) {
			allow_extend = 0;
		}


		// INDELS
		if(allow_diff && allow_indels) {
			if(e->state == STATE_I) {
				if(allow_extend) {
					// extend an insertion
					heap_push(heap, e->i-1, e->L, e->U, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps,  e->aln_length, e->aln_path, params);
				}
			} else {
				if(allow_open && (e->state == STATE_M)) {
					// open an insertion
					heap_push(heap, e->i-1, e->L, e->U, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
				for (int j = 1; j < alphabet_size; j++) {
					if (L[j] <= U[j]) {
						if(e->state == STATE_M) {
							if(allow_open) {
								// open a deletion
								heap_push(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo + 1, e->num_gape, STATE_D,
										e->num_snps, e->aln_length, e->aln_path, params);
							}
						} else {
							if(allow_extend) {
								// extend a deletion
								heap_push(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape + 1, STATE_D,
										e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
					}
				}
			}
		}

		// MATCH/MISMATCH
		int c = read[e->i-1];
		if(allow_diff && allow_mm) {
			for (int j = 1; j < alphabet_size; j++) {
				if (L[j] <= U[j]) { // check if this is a match or mismatch
					int is_mm = 0;
					if(is_multiref) {
						if((c > 3) || (j == 10)/*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0)) {
							is_mm = 1;
						}
					} else {
						if((c > 3) || (c != (j-1))) {
							is_mm = 1;
						}
					}
					if(!is_mm) {
						heap_push(heap, e->i-1, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					} else {
						heap_push(heap, e->i-1, L[j], U[j], e->num_mm+1, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
				}
			}
		} else if (c < 4) { // exact match only
			if(is_multiref) {
				for(int b = 0; b < BASES_PER_NUCLEOTIDE; b++) {
					unsigned char base = nucl_bases_table[c][b]; // does not match with N's in the reference
					if (L[base] <= U[base]) {
						heap_push(heap, e->i-1, L[base], U[base], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
				}
			} else {
				if (L[c+1] <= U[c+1]) {
					heap_push(heap, e->i-1, L[c+1], U[c+1], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
							e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
		}
	}
}

/* Priority Heap Operations */

priority_heap_t *heap_init(const aln_params_t *p) {
	priority_heap_t *heap = (priority_heap_t*) calloc(1, sizeof(priority_heap_t));
	// max number of possible different scores
	heap->num_buckets = aln_score(p->max_diff+1, p->max_gapo+1, p->max_gape+1, p);
	heap->buckets = (heap_bucket_t*) calloc(heap->num_buckets, sizeof(heap_bucket_t));
	if(heap == NULL || heap->buckets == NULL) {
		printf("Could not allocate memory for the heap \n");
		exit(1);
	}
	for (int i = 0; i < heap->num_buckets; i++) {
		heap_bucket_t *hb = &(heap->buckets[i]);
		hb->max_entries = 4;
		hb->entries = (aln_entry_t*) calloc(hb->max_entries, sizeof(aln_entry_t));
		if(hb->entries == NULL) {
			printf("Could not allocate memory for the heap \n");
			exit(1);
		}
	}
	heap->best_score = heap->num_buckets;
	return heap;
}

void heap_free(priority_heap_t *heap) {
	for (int i = 0; i < heap->num_buckets; i++) {
		free(heap->buckets[i].entries);
	}
	free(heap->buckets);
	free(heap);
}

void heap_reset(priority_heap_t *heap) {
	for (int i = 0; i < heap->num_buckets; i++) {
		heap->buckets[i].num_entries = 0;
	}
	heap->best_score = heap->num_buckets;
	heap->num_entries = 0;
}

void heap_push(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U,
							const int num_mm, const int num_gapo, const int num_gape,
							const int state, const int num_snps, const int aln_length, const char* aln_path,
							const aln_params_t *params) {
	
	// compute the alignment score for this entry
	int score = aln_score(num_mm, num_gapo, num_gape, params);
	heap_bucket_t *hb = &(heap->buckets[score]);
	if (hb->num_entries == hb->max_entries) {
		hb->max_entries <<= 1;
		hb->entries = (aln_entry_t*)realloc(hb->entries, sizeof(aln_entry_t) * hb->max_entries);
		if(hb->entries == NULL) {
			printf("Could not reallocate memory for the heap bucket! \n");
			exit(1);
        }
		// set the newly allocated memory to zero
		//memset(&(hb->entries[hb->num_entries]), 0, sizeof(aln_entry_t) * (hb->max_entries - hb->num_entries));
	}
	// new entry
	aln_entry_t *p = &(hb->entries[hb->num_entries]);
	p->i = i;
	p->score = score;
	p->L = L;
	p->U = U;
	p->num_mm = num_mm;
	p->num_gapo = num_gapo;
	p->num_gape = num_gape;
	p->state = state;
	p->num_snps = num_snps;
	p->aln_length = 0;
	if(aln_path != NULL) {
		memset(&(p->aln_path), 0, ALN_PATH_ALLOC*sizeof(char));
		memcpy(&(p->aln_path), aln_path, aln_length*sizeof(char));
		p->aln_path[aln_length] = state;
		p->aln_length = aln_length + 1;
	}

	hb->num_entries++;
	heap->num_entries++;

	if (heap->best_score > score) {
		heap->best_score = score;
	}
}

// returns the entry with the best (lowest) score
void heap_pop(priority_heap_t* heap, aln_entry_t* e) {
	heap_bucket_t* hb = &(heap->buckets[heap->best_score]);
	aln_entry_t* et = &(hb->entries[hb->num_entries - 1]);
	hb->num_entries--;
	heap->num_entries--;

	if ((hb->num_entries == 0) && heap->num_entries) { // find the next best score
		int i;
		for (i = heap->best_score + 1; i < heap->num_buckets; i++) {
			if (heap->buckets[i].num_entries != 0) break;
		}
		heap->best_score = i;
	} else if (heap->num_entries == 0) {
		heap->best_score = heap->num_buckets;
	}
	memcpy(e, et, sizeof(aln_entry_t));
}
