//-------------------------------------------------------
// Common Alignment Functionality
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
#include "exact_match.h"
#include "inexact_match.h"
#include "io.h"

void precalc_sa_intervals(bwt_t *BWT, const aln_params_t* params, char* preFname);
sa_intv_list_t* load_precalc_sa_intervals(const char* preFname);

void set_default_aln_params(aln_params_t* params) {
	params->gape_score = 4;
	params->gapo_score = 11;
	params->mm_score = 3;
	params->max_diff = 0;
	params->max_gape = 6;
	params->max_gapo = 1;
	params->seed_length = 32;
	params->max_diff_seed = 2;
	params->max_entries = 3000000;
	params->use_precalc = 0;
	params->matched_Ncontig = 0;
	params->is_multiref = 1;
	params->max_best = 30;
	params->no_indel_length = 5;
	params->n_threads = 1;
}

int align_reads(char* fastaFname, char* readsFname, char* alnsFname, aln_params_t* params) {
	printf("**** BWBBLE Read Alignment ****\n");
	char* bwtFname  = (char*) malloc(strlen(fastaFname) + 5);
	//char* alnsFname  = (char*) malloc(strlen(fastaFname) + 5);
	char* preFname = (char*) malloc(strlen(fastaFname) + 5);
	sprintf(bwtFname, "%s.bwt", fastaFname);
	//sprintf(alnsFname, "%s.aln", fastaFname);
	sprintf(preFname, "%s.pre", fastaFname);
	remove(alnsFname); // remove an older .aln file (if it exists)

	clock_t t = clock();
	bwt_t* BWT = load_bwt(bwtFname, 0);
	printf("Total BWT loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	t = clock();
	reads_t* reads = fastq2reads(readsFname);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);	

	sa_intv_list_t* sa_intv_table = NULL;
	if(params->use_precalc) {
		t = clock();
		if((FILE*) fopen(preFname, "r") == NULL) {
			precalc_sa_intervals(BWT, params, preFname);
		}
		sa_intv_table = load_precalc_sa_intervals(preFname);
		printf("Total pre-calculated intervals loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	t = clock();
	//if(params->max_diff == 0) {
		//align_reads_exact(BWT, reads, sa_intv_table, params, alnsFname);
	//} else {
		if(params->n_threads > 1) {
			align_reads_inexact_parallel(BWT, reads, sa_intv_table, params, alnsFname);
		} else {
			align_reads_inexact(BWT, reads, sa_intv_table, params, alnsFname);
		}
	//}
	printf("Total read alignment time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	free_bwt(BWT);
	free_reads(reads);
	if(sa_intv_table) free_sa_interval_list(sa_intv_table);
	free(bwtFname);
	//free(alnsFname);
	free(preFname);
	return 0;
}

/* SA Interval Management */

// intervals are always added in sorted order
// intervals cannot overlap but can be adjoining - adjoining intervals will be merged
void add_sa_interval(sa_intv_list_t* intv_list, bwtint_t L, bwtint_t U) {
	// adjoining intervals are merged
	if((intv_list->size != 0) && (L == (intv_list->last_intv->U + 1))) {
		intv_list->last_intv->U = U;
	} else {
		sa_intv_t* intv = (sa_intv_t*) calloc(1, sizeof(sa_intv_t));
		intv->L = L;
		intv->U = U;
		if(intv_list->size == 0) {
			intv_list->first_intv = intv;
			intv_list->last_intv = intv;
		} else {
			intv_list->last_intv->next_intv = intv;
			intv_list->last_intv = intv;
		}
		intv_list->size++;
	}
}

void clear_sa_interval_list(sa_intv_list_t* intv_list) {
	sa_intv_t* intv = intv_list->first_intv;
	for(int s = 0; s < intv_list->size; s++) {
		sa_intv_t* tmp = intv;
		intv = intv->next_intv;
		free(tmp);
	}
	intv_list->size = 0;
	intv_list->first_intv = 0;
	intv_list->last_intv = 0;
}

void free_sa_interval_list(sa_intv_list_t* intv_list) {
	sa_intv_t* intv = intv_list->first_intv;
	for(int s = 0; s < intv_list->size; s++) {
		sa_intv_t* tmp = intv;
		intv = intv->next_intv;
		free(tmp);
	}
	free(intv_list);
}

void print_sa_interval_list(sa_intv_list_t* intv_list) {
	if (intv_list->size != 0) {
		sa_intv_t* intv = intv_list->first_intv;
		for(int i = 0; i < intv_list->size; i++) {
			printf("SA Interval %d: L = %" PRIbwtint_t ", U = %" PRIbwtint_t "\n", i, intv->L, intv->U);
			intv = intv->next_intv;
		}
	}
}

void store_sa_interval_list(sa_intv_list_t* intv_list, FILE* saFile) {
	fwrite(&intv_list->size, sizeof(int), 1, saFile);
	sa_intv_t* intv = intv_list->first_intv;
	for(int s = 0; s < intv_list->size; s++) {
		fwrite(&intv->L, sizeof(bwtint_t), 1, saFile);
		fwrite(&intv->U, sizeof(bwtint_t), 1, saFile);
		intv = intv->next_intv;
	}
}

void load_sa_intervals_error() {
	printf("load_sa_intervals: Could not read the precomputed intervals from file! \n");
	exit(1);
}

void load_sa_interval_list(sa_intv_list_t* intv_list, FILE* saFile) {
	if(fread(&intv_list->size, sizeof(int), 1, saFile) < 1) load_sa_intervals_error();
	for(int j = 0; j < intv_list->size; j++) {
		sa_intv_t* intv = (sa_intv_t*) calloc(1, sizeof(sa_intv_t));
		if(fread(&intv->L, sizeof(bwtint_t), 1, saFile) < 1) load_sa_intervals_error();
		if(fread(&intv->U, sizeof(bwtint_t), 1, saFile) < 1) load_sa_intervals_error();
		if(j == 0) {
			intv_list->first_intv = intv;
			intv_list->last_intv = intv;
		} else {
			intv_list->last_intv->next_intv = intv;
			intv_list->last_intv = intv;
		}
	}
}
int read2index(char* read, int readLen) {
  int index = 0;
  for(int i = readLen-PRECALC_INTERVAL_LENGTH; i < readLen; i++) {
	  if(read[i] >= NUM_NUCLEOTIDES) {
		  // N's are treated as mismatches
		  return -1;
	  }
	  index *= NUM_NUCLEOTIDES;
	  index += read[i];
  }
  return index;
}
void next_read(read_t* read) {
	read->seq[read->len-1]++;
	for(int i = read->len - 1; i > 0; i--) {
		if(read->seq[i] < NUM_NUCLEOTIDES) {
			break;
		}
		read->seq[i] -= NUM_NUCLEOTIDES;
		read->seq[i-1]++;
	}
	if(read->seq[0] >= NUM_NUCLEOTIDES) {
		read->seq[0] -= NUM_NUCLEOTIDES;
	}
}

void precalc_sa_intervals(bwt_t *BWT, const aln_params_t* params, char* preFname) {
	printf("Pre-calculating SA intervals...\n");
	FILE* preFile = (FILE*)fopen(preFname, "wb");
	if (preFile == NULL) {
		fprintf(stderr, "precalc_sa_intervals: Cannot open PRE file %s!\n", preFname);
		exit(1);
	}
	sa_intv_list_t** sa_intervals = (sa_intv_list_t**) malloc(NUM_PRECALC*sizeof(sa_intv_list_t*));
	read_t *read = (read_t*) calloc(1, sizeof(read_t));
	read->len = PRECALC_INTERVAL_LENGTH;
	read->seq = (char*) calloc(read->len, sizeof(unsigned int));
	clock_t t = clock();
	for(int i = 0; i < NUM_PRECALC; i++) {
		exact_match(BWT, read, &sa_intervals[i], params);
		next_read(read);
	}
	printf("Interval pre-computation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	t = clock();
	for(int i = 0; i < NUM_PRECALC; i++) {
		store_sa_interval_list(sa_intervals[i], preFile);
		free_sa_interval_list(sa_intervals[i]);
	}
	fclose(preFile);
	printf("Storing results time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}

sa_intv_list_t* load_precalc_sa_intervals(const char* preFname) {
	FILE* preFile = (FILE*)fopen(preFname, "rb");
	if (preFile == NULL) {
		fprintf(stderr, "load_sa_intervals: Cannot open the PRE file: %s!\n", preFname);
		exit(1);
	}
	sa_intv_list_t* intv_list = (sa_intv_list_t*) malloc(NUM_PRECALC * sizeof(sa_intv_list_t));
	for(int i = 0; i < NUM_PRECALC; i++) {
		load_sa_interval_list(&intv_list[i], preFile);
	}
	fclose(preFile);
	return intv_list;
}

/* Alignments Management */

static inline int aln_score(const int m, const int o, const int e, const aln_params_t *p) {
	return m*p->mm_score + o*p->gapo_score + e*p->gape_score;
}

// allocate memory to store the alignments
alns_t* init_alignments() {
	alns_t* alns = (alns_t*) calloc(1, sizeof(alns_t));
	alns->num_entries = 0;
	alns->max_entries = 4;
	alns->entries = (aln_t*) calloc(alns->max_entries, sizeof(aln_t));
	if(alns == NULL || alns->entries == NULL) {
		printf("Could not allocate memory for alignments \n");
		exit(1);
	}
	return alns;
}

void free_alignments(alns_t* alns) {
	for(int i = 0; i < alns->num_entries; i++) {
		if(alns->entries[i].aln_path) free(alns->entries[i].aln_path);
	}
	free(alns->entries);
	free(alns);
}

void reset_alignments(alns_t* alns) {
	alns->num_entries = 0;
}

void add_alignment(aln_entry_t* e, const bwtint_t L, const bwtint_t U, int score, alns_t* alns, const aln_params_t* params) {
	// do not add the alignment if we already have an alignment with these bounds (can occur with gaps)
	if(e->num_gapo) {
		for(int j = 0; j < alns->num_entries; j++) {
			aln_t aln = alns->entries[j];
			if(aln.L == L && aln.U == U) {
				return;
			}
		}
	}
	if(alns->num_entries == alns->max_entries) {
		alns->max_entries <<= 1;
		alns->entries = (aln_t*)realloc(alns->entries, alns->max_entries * sizeof(aln_t));
		memset(alns->entries + alns->max_entries/2, 0,  (alns->max_entries/2)*sizeof(aln_t));
	}
	aln_t *alt = &(alns->entries[alns->num_entries]);
	alt->num_mm = e->num_mm;
	alt->num_gapo = e->num_gapo;
	alt->num_gape = e->num_gape;
	alt->num_snps = e->num_snps;
	alt->L = L;
	alt->U = U;
	alt->score = score;
	alt->aln_length = e->aln_length;
	alt->aln_path = (char*) malloc(e->aln_length*sizeof(char));
	memcpy(alt->aln_path, &(e->aln_path), e->aln_length*sizeof(char));
	alns->num_entries++;
}

void print_alignments(alns_t* alns) {
	printf("Number of alignments = %d \n", alns->num_entries);
	for(int j = 0; j < alns->num_entries; j++) {
		aln_t aln = alns->entries[j];
		printf("Alignment %d: SA(%" PRIbwtint_t ",%" PRIbwtint_t ") score = %d, num_mm = %u, num_go = %u, num_ge = %u, num_snps = %u, aln_length = %d\n",
				j, aln.L, aln.U, aln.score, aln.num_mm, aln.num_gapo, aln.num_gape, aln.num_snps, 0);//, aln.aln_length);
		printf("\n");
	}
}

// create alignments for SA intervals from exact matching
alns_t* sa_intervals2alns(sa_intv_list_t* intv_list, int aln_length) {
	alns_t* alns = (alns_t*) calloc(1, sizeof(alns_t));
	alns->num_entries = intv_list->size;
	alns->entries = (aln_t*) calloc(alns->num_entries, sizeof(aln_t));
	if(alns == NULL || alns->entries == NULL) {
		printf("Could not allocate memory for alignments\n");
		exit(1);
	}
	sa_intv_t* intv = intv_list->first_intv;
	for(int i = 0; i < intv_list->size; i++) {
		aln_t* aln = &(alns->entries[i]);
		aln->L = intv->L;
		aln->U = intv->U;
		aln->aln_length = aln_length;
		// all other aln parameters are 0
		intv = intv->next_intv;
	}
	return alns;
}

// store alignments to file
void alns2alnf(alns_t* alns, FILE* alnFile) {
	fprintf(alnFile, "%d\n", alns->num_entries);
	for(int i = 0; i < alns->num_entries; i++) {
		aln_t aln = alns->entries[i];
		fprintf(alnFile, "%d\t%llu\t%llu\t%d\t%d\t%d\t%d\t", aln.score, (unsigned long long int) aln.L, (unsigned long long int) aln.U,
				aln.num_mm, aln.num_gapo, aln.num_gape, /*aln.num_snps,*/ aln.aln_length);
		for(int j = aln.aln_length-1; j >= 0; j--) {
			fprintf(alnFile, "%c ", aln.aln_path[j]);
		}
		fprintf(alnFile, "\n");
	}
}

void alns2alnf_bin(alns_t* alns, FILE* alnFile) {
	fwrite(&alns->num_entries, sizeof(int), 1, alnFile);
	for(int i = 0; i < alns->num_entries; i++) {
		aln_t aln = alns->entries[i];
		fwrite(&aln.score, sizeof(int), 1, alnFile);
		fwrite(&aln.L, sizeof(unsigned long long int), 1, alnFile);
		fwrite(&aln.U, sizeof(unsigned long long int), 1, alnFile);
		fwrite(&aln.num_mm, sizeof(int), 1, alnFile);
		fwrite(&aln.num_gapo, sizeof(int), 1, alnFile);
		fwrite(&aln.num_gape, sizeof(int), 1, alnFile);
		fwrite(&aln.aln_length, sizeof(int), 1, alnFile);

		int state_pairs = 0;
		if(aln.aln_length > 0) {
			int* states = (int*) calloc(aln.aln_length, sizeof(int));
			int state = aln.aln_path[aln.aln_length-1];
			uint16_t state_counter = 1;
			state_pairs = 1;
			for(int j = aln.aln_length-2; j >= 0; j--) {
				if(state == aln.aln_path[j]) {
					state_counter++;
				} else {
					states[state_pairs - 1] = state | (state_counter << 2);
					state = aln.aln_path[j];
					state_counter = 1;
					state_pairs++;
				}
			}
			states[state_pairs-1] = state | (state_counter << 2);
			fwrite(&state_pairs, sizeof(int), 1, alnFile);
			for(int j = 0; j < state_pairs; j++) {
				fwrite(&states[j], sizeof(int), 1, alnFile);
			}
		} else {
			fwrite(&state_pairs, sizeof(int), 1, alnFile);
		}
	}
}

// load alignments from file

void load_alns_error(const char* alnFname) {
	printf("alnsf2alns: Could not parse read alignment data in file (file content did not match expected format): %s!\n", alnFname);
	exit(1);
}

alns_t* alnsf2alns(int* n_alns, char *alnFname) {
	FILE * alnFile = (FILE*) fopen(alnFname, "r");
	if (alnFile == NULL) {
		printf("alnsf2alns: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}
	int num_alns = 0;
	int alloc_alns = 2000000;
	alns_t* alns = (alns_t*) calloc(alloc_alns, sizeof(alns_t));

	while(!feof(alnFile)) {
		if(num_alns == alloc_alns) {
			alloc_alns <<= 1;
			alns = (alns_t*) realloc(alns, alloc_alns*sizeof(alns_t));
			memset(alns + alloc_alns/2, 0,  (alloc_alns/2)*sizeof(alns_t));
		}
		alns_t* read_alns = &alns[num_alns];
		if(fscanf(alnFile, "%d\n", &(read_alns->num_entries)) != 1) load_alns_error(alnFname);
		read_alns->entries = (aln_t*) calloc(read_alns->num_entries, sizeof(aln_t));
		for(int i = 0; i < read_alns->num_entries; i++) {
			aln_t* aln = &(read_alns->entries[i]);
			if(fscanf(alnFile, "%d\t%" SCNbwtint_t "\t%" SCNbwtint_t "\t%d\t%d\t%d\t%d\t", &(aln->score), &(aln->L), &(aln->U),
					&(aln->num_mm), &(aln->num_gapo), &(aln->num_gape), /*&(aln->num_snps),*/ &(aln->aln_length)) != 7) {
				load_alns_error(alnFname);
			}
			aln->aln_path = (char*) malloc(aln->aln_length*sizeof(char));
			for(int j = 0; j < aln->aln_length; j++) {
				if(fscanf(alnFile, "%c ", &(aln->aln_path[j])) != 1) load_alns_error(alnFname);
			}
			if(fscanf(alnFile,"\n") < 0) load_alns_error(alnFname);
		}
		num_alns++;
	}
	*n_alns = num_alns;
	fclose(alnFile);
	return alns;
}

alns_t* alnsf2alns_bin(int* n_alns, char *alnFname) {
	FILE * alnFile = (FILE*) fopen(alnFname, "rb");
	if (alnFile == NULL) {
		printf("alnsf2alns: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}
	int num_alns = 0;
	int alloc_alns = 2000000;
	alns_t* alns = (alns_t*) calloc(alloc_alns, sizeof(alns_t));

	while(!feof(alnFile)) {
		if(num_alns == alloc_alns) {
			alloc_alns <<= 1;
			alns = (alns_t*) realloc(alns, alloc_alns*sizeof(alns_t));
			memset(alns + alloc_alns/2, 0,  (alloc_alns/2)*sizeof(alns_t));
		}
		alns_t* read_alns = &alns[num_alns];
		if(fread(&(read_alns->num_entries), sizeof(int), 1, alnFile) < 1) {
			if(feof(alnFile)) break;
			load_alns_error(alnFname);
		}
		read_alns->entries = (aln_t*) calloc(read_alns->num_entries, sizeof(aln_t));
		for(int i = 0; i < read_alns->num_entries; i++) {
			aln_t* aln = &(read_alns->entries[i]);

			if(fread(&(aln->score), sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
			if(fread(&(aln->L), sizeof(unsigned long long int), 1, alnFile) < 1) load_alns_error(alnFname);
			if(fread(&(aln->U), sizeof(unsigned long long int), 1, alnFile) < 1) load_alns_error(alnFname);
			if(fread(&(aln->num_mm), sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
			if(fread(&(aln->num_gapo), sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
			if(fread(&(aln->num_gape), sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
			if(fread(&(aln->aln_length), sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
			aln->aln_path = (char*) malloc(aln->aln_length*sizeof(char));
			int state_pairs;
			if(fread(&state_pairs, sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
			int state_pair;
			int path_idx = 0;
			for(int j = 0; j < state_pairs; j++) {
				if(fread(&state_pair, sizeof(int), 1, alnFile) < 1) load_alns_error(alnFname);
				int counter = state_pair >> 2;
				int state = state_pair & 3; // last 2 bits
				for(int k = 0; k < counter; k++) {
					aln->aln_path[path_idx] = state;
					path_idx++;
				}
			}
		}
		num_alns++;
	}
	*n_alns = num_alns;
	fclose(alnFile);
	return alns;
}

/* Alignment Result Evaluation & SAM IO */

// SAM IO / MAPQ adapted from BWA

int check_ref_mapping(read_t* read, int is_multiref);
void eval_aln(read_t* read, alns_t* alns, bwt_t* BWT, int is_multiref, int max_mm);
void print_aln2sam(FILE* samFile, const fasta_annotations_t* annotations, read_t* r);


void alns2sam(char *fastaFname, char *readsFname, char *alnsFname, char* samFname, int is_multiref, int max_diff) {
	printf("**** BWBBLE Alignment Evaluation/SAM File Generation ****\n");

	// load the BWT
	char* bwtFname  = (char*) malloc(strlen(fastaFname) + 5);
	char* annFname  = (char*) malloc(strlen(fastaFname) + 5);
	sprintf(bwtFname, "%s.bwt", fastaFname);
	sprintf(annFname, "%s.ann", fastaFname);
	bwt_t* BWT = load_bwt(bwtFname, 1);
	fasta_annotations_t* annotations = annf2ann(annFname);

	// load the alignment results of all the reads (TODO: batch)
	int num_alns;
	alns_t* alns = alnsf2alns_bin(&num_alns, alnsFname);
	//alns_t* alns = alnsf2alns(&num_alns, alnsFname);
	reads_t* reads = fastq2reads(readsFname);
	//assert(num_alns == reads->count);

	// open SAM for writing
	FILE* samFile = (FILE*) fopen(samFname, "w");
	if (samFile == NULL) {
		printf("alns2sam: Cannot open SAM file: %s!\n", samFname);
		perror(samFname);
		exit(1);
	}

	// print SAM headers
	// @SQ (name, length)
	for (int i = 0; i < annotations->num_seq; i++) {
		// TODO: avoid printing bubble branches
		fprintf(samFile, "@SQ\tSN:%s\tLN:%d\n", annotations->seq_anns[i].name, (int) (annotations->seq_anns[i].end_index - annotations->seq_anns[i].start_index+1));
	}
	// @RG
	// @PG
	fprintf(samFile, "@PG\tID:bwbble\tPN:bwbble\tVN:0.1-r01\n");

	// evaluate alignment results per read / output in SAM format
	int num_processed = 0;
	while(num_processed < reads->count) {
		int batch_size = ((reads->count - num_processed) > READ_BATCH_SIZE ) ? READ_BATCH_SIZE : (reads->count - num_processed);
		for(int i = 0; i < batch_size; i++) {
			if(num_processed == num_alns) {
				break;
			}
			read_t* read = &reads->reads[num_processed];
			eval_aln(read, &alns[num_processed], BWT, is_multiref, max_diff);
			print_aln2sam(samFile, annotations, read);
			num_processed++;
		}
		printf("Processed %d reads.\n", num_processed);
		if(num_processed == num_alns) {
			break;
		}
	}

	free(bwtFname);
	free(annFname);
	free_bwt(BWT);
	free_reads(reads);
	free_alignments(alns);
	free_ann(annotations);
	fclose(samFile);
}

// BWA Flags
#define SAM_FSU   4 // self-unmapped
#define SAM_FSR  16 // self on the reverse strand

void print_aln2sam(FILE* samFile, const fasta_annotations_t* annotations, read_t* r) {
	int flag = 0; // FLAG
	if(r->aln_type != ALN_NOMATCH) {
		int seqid = -1;
		for(int i = 0; i < annotations->num_seq; i++) { // TODO: binary search
			if((r->aln_pos >= annotations->seq_anns[i].start_index) && (r->aln_pos <= annotations->seq_anns[i].end_index)) {
				seqid = i; break;
			}
		}
		if (r->aln_strand) flag |= SAM_FSR;

		// QNAME, FLAG, RNAME
		fprintf(samFile, "%s\t%d\t%s\t", r->name, flag, annotations->seq_anns[seqid].name);
		// POS (1-based), MAPQ
		fprintf(samFile, "%d\t%d\t", (int)(r->aln_pos - annotations->seq_anns[seqid].start_index + 1), r->mapQ);

		// CIGAR
		if (r->aln_strand) { // reverse aln path
			for (int i = 0; i < r->aln_length >> 1; i++) {
				char tmp = r->aln_path[r->aln_length-1-i];
				r->aln_path[r->aln_length-1-i] = r->aln_path[i]; r->aln_path[i] = tmp;
			}
		}
		int cigar_length = 1;
		unsigned char last_type = r->aln_path[0];
		for (int i = 1; i < r->aln_length; i++) {
			if (last_type != r->aln_path[i]) cigar_length++;
			last_type = r->aln_path[i];
		}
		int* cigar = (int*)malloc(cigar_length * sizeof(int));
		cigar[0] = 1u << 4 | r->aln_path[r->aln_length-1];
		last_type = r->aln_path[r->aln_length-1];

		int cigar_idx = 0;
		for (int i = r->aln_length - 2; i >= 0; i--) {
			if (r->aln_path[i] == last_type) {
				cigar[cigar_idx] += 1u << 4;
			} else {
				cigar_idx++;
				cigar[cigar_idx] = 1u << 4 | r->aln_path[i];
				last_type = r->aln_path[i];
			}
		}
		for (int i = 0; i < cigar_length; i++) {
			fprintf(samFile, "%d%c", (cigar[i] >> 4), "MID"[(cigar[i] & 0xf)]);
		}

		// RNEXT, PNEXT, TLEN (print void mate position and coordinate)
		fprintf(samFile, "\t*\t0\t0\t");

		// SEQ, QUAL (print sequence and quality)
		char* seq = r->aln_strand ? r->rc : r->seq;
		for (int i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", "AGCTN"[(int)seq[i]]);
		}
		fprintf(samFile, "\t");

		if (r->qual) {
			if (r->aln_strand) { // reverse quality
				for (int i = 0; i < r->len >> 1; i++) {
					char tmp = r->qual[r->len-1-i];
					r->qual[r->len-1-i] = r->qual[i]; r->qual[i] = tmp;
				}
			}
			fprintf(samFile, "%s", r->qual);
		} else fprintf(samFile, "*");
		fprintf(samFile, "\n");

	} else { // unmapped read
		char* seq = r->aln_strand ? r->rc : r->seq;
		int flag = SAM_FSU;
		// QNAME, FLAG
		fprintf(samFile, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", r->name, flag);

		// SEQ, QUAL (print sequence and quality)
		for (int i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", "AGCTN"[(int)seq[i]]);
		}
		fprintf(samFile, "\t");
		if (r->qual) {
			if (r->aln_strand) { // reverse quality
				for (int i = 0; i < r->len >> 1; i++) {
					char tmp = r->qual[r->len-1-i];
					r->qual[r->len-1-i] = r->qual[i]; r->qual[i] = tmp;
				}
			}
			fprintf(samFile, "%s", r->qual);
		} else fprintf(samFile, "*");
		fprintf(samFile, "\n");
	}
}

// Evaluate the alignment results of a given set of reads
void eval_alns(char *fastaFname, char *readsFname, char *alnsFname, int is_multiref, int max_diff) {
	printf("**** BWBBLE Alignment Evaluation ****\n");

	// categorize reads (ids) by alignment type
	FILE* unalignedFile = (FILE*) fopen("bwbble.unaligned", "wb");
	FILE* confidentFile = (FILE*) fopen("bwbble.conf", "wb");
	FILE* correctFile = (FILE*) fopen("bwbble.corr", "wb");
	FILE* misalignedFile = (FILE*) fopen("bwbble.mis", "wb");
	if ((unalignedFile == NULL) || (confidentFile == NULL) || (correctFile == NULL) || (misalignedFile == NULL)) {
		printf("eval: Cannot open the unaligned/conf/corr/mis file(s)!\n");
		exit(1);
	}

	char* bwtFname  = (char*) malloc(strlen(fastaFname) + 5);
	sprintf(bwtFname, "%s.bwt", fastaFname);

	// load the alignment results of all the reads
	int num_alns;
	alns_t* alns = alnsf2alns(&num_alns, alnsFname);
	bwt_t* BWT = load_bwt(bwtFname, 1);
	reads_t* reads = fastq2reads(readsFname);
	assert(num_alns == reads->count);

	// for each read: evaluate the alignment quality and accuracy
	int n_confident = 0;
	int n_correct = 0;
	int n_misaligned = 0;
	int n_unaligned = 0;

	for(int i = 0; i < reads->count; i++) {
		read_t* read = &(reads->reads[i]);
		parse_read_mapping(read);
		eval_aln(read, &alns[i], BWT, is_multiref, max_diff);
		if(read->aln_type == ALN_NOMATCH) {
			n_unaligned++;
			fwrite(&i, sizeof(int), 1, unalignedFile);
			continue;
		}
		if(read->mapQ < MAPQ_CONFIDENT) {
			continue;
		}
		n_confident++;
		fwrite(&i, sizeof(int), 1, confidentFile);
		if(check_ref_mapping(read, is_multiref) == 1) {
			n_correct++;
			fwrite(&i, sizeof(int), 1, correctFile);
		} else {
			n_misaligned++;
			fwrite(&i, sizeof(int), 1, misalignedFile);
		}
	}
	fwrite(&n_unaligned, sizeof(int), 1, unalignedFile);
	fwrite(&n_confident, sizeof(int), 1, confidentFile);
	fwrite(&n_correct, sizeof(int), 1, correctFile);
	fwrite(&n_misaligned, sizeof(int), 1, misalignedFile);

	printf("total num_reads = %d, confident = %d correct = %d, misaligned = %d, unaligned = %d\n", reads->count, n_confident, n_correct, n_misaligned, n_unaligned);

	free(bwtFname);
	free_reads(reads);
	free_alignments(alns);
	free_bwt(BWT);

	fclose(unalignedFile);
	fclose(confidentFile);
	fclose(correctFile);
	fclose(misalignedFile);
}

int mapq2(const read_t *read, int max_mm, int is_multiref) {
	if (read->aln_top1_count == 0) return 23; // no hits
	if(is_multiref) {
		if (read->aln_top1_count > read->num_mref_pos) return 0; // repetitive top hit
	} else {
		if (read->aln_top1_count > (read->ref_pos_r - read->ref_pos_l + 1)) return 0;
	}
	if (read->num_mm == max_mm) return 25;
	if (read->aln_top2_count == 0) return 37; // unique
	int n = (read->aln_top2_count >= 255)? 255 : read->aln_top2_count;
	int q = (int)(4.343 * log(n) + 0.5);
	return (23 < q)? 0 : 23 - q;
}

int mapq(const read_t *read, int max_mm, int is_multiref) {
	if (read->aln_top1_count == 0) return 23; // no hits
	if (read->aln_top1_count > 1) return 0; // repetitive top hit
	if (read->num_mm == max_mm) return 25;
	if (read->aln_top2_count == 0) return 37; // unique
	int n = (read->aln_top2_count >= 255)? 255 : read->aln_top2_count;
	int q = (int)(4.343 * log(n) + 0.5);
	return (23 < q)? 0 : 23 - q;
}

int get_aln_length(char* aln_path, int path_length) {
	int aln_length = path_length;
	// discard all the insertions
	for(int i = 0; i < path_length; i++) {
		if(aln_path[i] == STATE_I) {
			aln_length--;
		}
	}
	return aln_length;
}

// Evaluate the alignment results of a given read
void eval_aln(read_t* read, alns_t* alns, bwt_t* BWT, int is_multiref, int max_mm) {
	// no matches
	if(alns->num_entries == 0) {
		read->aln_top1_count = 0;
		read->aln_top2_count = 0;
		read->aln_type = ALN_NOMATCH;
		return;
	}

	int best_score = alns->entries[0].score;
	for(int i = 0; i < alns->num_entries; i++) {
		aln_t aln = alns->entries[i];
		if(aln.score > best_score) {
			read->aln_top2_count += (aln.U - aln.L + 1);
		} else {
			// select one of the top score alignments
			read->aln_top1_count += (aln.U - aln.L + 1);
			// pick only 1 top aln for this read
			if(i == 0) {
				read->num_mm = aln.num_mm;
				read->num_gapo = aln.num_gapo;
				read->num_gape = aln.num_gape;
				read->aln_score = aln.score;
				read->aln_length = aln.aln_length;
				read->aln_path = (char*) malloc(aln.aln_length*sizeof(char));
				memcpy(read->aln_path, aln.aln_path, aln.aln_length*sizeof(char));
				// pick one of the matches (TODO: pick randomly from the L,U range)
				read->aln_sa = aln.L;// + (bwtint_t)((aln.U - aln.L + 1));
				// determine the position and strand of the mapping
				bwtint_t ref_pos = SA(BWT, read->aln_sa);
				if(ref_pos > (BWT->length-1)/2) {
					read->aln_strand = 0; // read rc + ref rc <=> read fwd + ref fwd
					bwtint_t fwd_pos = (BWT->length - 1) - ref_pos - 1;
					read->aln_pos = fwd_pos - get_aln_length(aln.aln_path, aln.aln_length) + 1;
				} else {
					read->aln_strand = 1; // read rc + fwd ref <=> fwd read/ref rc
					//bwtint_t rc_pos = (BWT->length - 1) - ref_pos - 1;
					//read->aln_pos = rc_pos - get_aln_length(aln.aln_path, aln.aln_length) + 1 - (BWT->length-1)/2;
					read->aln_pos = ref_pos;
				}
			}
		}
	}

	/*if(is_multiref) {
		read->aln_type = (read->aln_top1_count > read->num_mref_pos) ? ALN_REPEAT : ALN_UNIQUE;
	} else {
		read->aln_type = (read->aln_top1_count > read->ref_pos_r - read->ref_pos_l + 1) ? ALN_REPEAT : ALN_UNIQUE;
	}*/

	read->aln_type = (read->aln_top1_count > 1) ? ALN_REPEAT : ALN_UNIQUE;
	read->mapQ = mapq(read, max_mm, is_multiref);
}

// 1 - correct, 0 - incorrect
int check_ref_mapping(read_t* read, int is_multiref) {
	// check strands
	if((read->aln_strand && !read->strand) || (!read->aln_strand && read->strand)) {
		return 0;
	}
	// check position
	if(is_multiref) {
		for(int i = 0; i < read->num_mref_pos; i++) {
			if(read->aln_pos == read->mref_pos[i] - 1) {
				return 1;
			}
		}
	} else {
		for(bwtint_t i = 0; i < read->ref_pos_r - read->ref_pos_l + 1; i++) {
			if(read->aln_pos == (read->ref_pos_l + i) - 1) {
				return 1;
			}
		}
	}
	return 0;
}
