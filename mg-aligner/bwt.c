//------------------------------------------------------
// BWT Indexing
// Victoria Popic (viq@stanford.edu), 2012
//-------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "bwt.h"
#include "io.h"

void compute_C(bwt_t *BWT, const unsigned char *seqIUPAC);
void compute_O(bwt_t *BWT);
void compute_SA(bwt_t *BWT);
void generate_occ_table(bwt_t *BWT);
bwtint_t get_occ_count(const bwt_t *BWT, const unsigned char c, const bwtint_t start_index, const bwtint_t end_index);
inline bwtint_t get_occ_count_opt(const bwt_t *BWT, const unsigned char c, const bwtint_t start_index, const bwtint_t end_index);
bwtint_t get_occ_count_opt_sse(const bwt_t *BWT, const unsigned char c, const bwtint_t start_index, const bwtint_t end_index);
inline void get_occ_count_alphabet(const bwt_t *BWT, const bwtint_t start_index, const bwtint_t end_index, int alphabet_size, bwtint_t *occ);
inline void get_occ_count_actg_alphabet(const bwt_t *BWT, const bwtint_t start_index, const bwtint_t end_index, bwtint_t *occ);

// Compute the BWT index of a given reference file
int index_bwt(const char *fastaFname, const char *extSAFname)
{
	printf("**** BWT Index **** \n");
	char *annFname = (char *)malloc(strlen(fastaFname) + 5);
	char *bwtFname = (char *)malloc(strlen(fastaFname) + 5);
	char *refFname = (char *)malloc(strlen(fastaFname) + 5);
	sprintf(annFname, "%s.ann", fastaFname);
	sprintf(bwtFname, "%s.bwt", fastaFname);
	sprintf(refFname, "%s.ref", fastaFname);

	// 1. read and pack the input FASTA file
	unsigned char *seq;
	bwtint_t seqLen;
	// resulting sequence will contain the forward + reverse complement reference
	if (extSAFname == NULL)
	{
		fasta2ref(fastaFname, refFname, annFname, &seq, &seqLen);
	}
	else
	{
		// load directly from precomputed ref file
		ref2seq(refFname, &seq, &seqLen);
	}

	// 2. compute the forward and reverse complement BWT
	printf("Running BWT construction\n");
	clock_t t = clock();
	bwt_t *BWT = construct_bwt(seq, seqLen, extSAFname);
	printf("Total BWT construction time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	//print_bwt(BWT);

	// 3. save the BWT to file
	store_bwt(BWT, bwtFname);

	free(annFname);
	free(bwtFname);
	free(refFname);
	free_bwt(BWT);
	return 0;
}

// Store the BWT to a file
void store_bwt(const const bwt_t *BWT, const char *bwtFname)
{
	FILE *bwtFile = (FILE *)fopen(bwtFname, "wb");
	if (bwtFile == NULL)
	{
		printf("store_bwt: Cannot open the BWT file %s!\n", bwtFname);
		fflush(stdout);
		exit(1);
	}
	fwrite(&BWT->length, sizeof(bwtint_t), 1, bwtFile);
	fwrite(&BWT->num_words, sizeof(bwtint_t), 1, bwtFile);
	fwrite(&BWT->num_sa, sizeof(bwtint_t), 1, bwtFile);
	fwrite(&BWT->num_occ, sizeof(bwtint_t), 1, bwtFile);
	fwrite(&BWT->sa0_index, sizeof(bwtint_t), 1, bwtFile);
	fwrite(&BWT->C, sizeof(bwtint_t), ALPHABET_SIZE + 1, bwtFile);
	fwrite(BWT->bwt, sizeof(uint32_t), BWT->num_words, bwtFile);
	fwrite(BWT->O, sizeof(bwtint_t), BWT->num_occ * ALPHABET_SIZE, bwtFile);
	fwrite(BWT->SA, sizeof(bwtint_t), BWT->num_sa, bwtFile);
	fclose(bwtFile);
}

void load_bwt_error(const char *bwtFname)
{
	printf("load_bwt: Could not read BWT from file: %s!\n", bwtFname);
	fflush(stdout);
	exit(1);
}

// Load the BWT from a file (for alignment only, SA not necessary)
bwt_t *load_bwt(const char *bwtFname, int loadSA)
{
	FILE *bwtFile = (FILE *)fopen(bwtFname, "rb");
	if (bwtFile == NULL)
	{
		printf("load_bwt: Cannot open the BWT file: %s!\n", bwtFname);
		fflush(stdout);
		exit(1);
	}
	bwt_t *BWT = (bwt_t *)calloc(1, sizeof(bwt_t));
	if (fread(&BWT->length, sizeof(bwtint_t), 1, bwtFile) < 1)
		load_bwt_error(bwtFname);
	if (fread(&BWT->num_words, sizeof(bwtint_t), 1, bwtFile) < 1)
		load_bwt_error(bwtFname);
	if (fread(&BWT->num_sa, sizeof(bwtint_t), 1, bwtFile) < 1)
		load_bwt_error(bwtFname);
	if (fread(&BWT->num_occ, sizeof(bwtint_t), 1, bwtFile) < 1)
		load_bwt_error(bwtFname);
	if (fread(&BWT->sa0_index, sizeof(bwtint_t), 1, bwtFile) < 1)
		load_bwt_error(bwtFname);
	if (fread(&BWT->C, sizeof(bwtint_t), ALPHABET_SIZE + 1, bwtFile) < ALPHABET_SIZE + 1)
		load_bwt_error(bwtFname);
	//BWT->bwt = (uint32_t *) _mm_malloc(BWT->num_words*sizeof(uint32_t), 16);
	BWT->bwt = (uint32_t *)calloc(BWT->num_words, sizeof(uint32_t));
	BWT->O = (bwtint_t *)calloc(BWT->num_occ * ALPHABET_SIZE, sizeof(bwtint_t));
	if ((BWT->bwt == 0) || (BWT->O == 0))
	{
		printf("load_bwt: Could not allocate memory for the BWT index. \n");
		fflush(stdout);
		exit(1);
	}
	if (fread(BWT->bwt, sizeof(uint32_t), BWT->num_words, bwtFile) < BWT->num_words)
		load_bwt_error(bwtFname);
	if (fread(BWT->O, sizeof(bwtint_t), BWT->num_occ * ALPHABET_SIZE, bwtFile) < BWT->num_occ * ALPHABET_SIZE)
		load_bwt_error(bwtFname);

	if (loadSA != 0)
	{
		BWT->SA = (bwtint_t *)calloc(BWT->num_sa, sizeof(bwtint_t));
		if (BWT->SA == 0)
		{
			printf("load_bwt: Could not allocate memory for the BWT index. \n");
			fflush(stdout);
			exit(1);
		}
		if (fread(BWT->SA, sizeof(bwtint_t), BWT->num_sa, bwtFile) < BWT->num_sa)
			load_bwt_error(bwtFname);
	}
	generate_occ_table(BWT);

	fclose(bwtFile);
	return BWT;
}

void load_ext_sa_error(const char *extSAFname)
{
	printf("esa2bwt: Could not read ext SA from file: %s!\n", extSAFname);
	fflush(stdout);
	exit(1);
}

bwtint_t esa2bwt(const unsigned char *seq, unsigned char *bwt_seq, const bwtint_t n, bwtint_t *bwtSA, const char *extSAFname)
{
	FILE *saFile = (FILE *)fopen(extSAFname, "rb");
	if (saFile == NULL)
	{
		printf("esa2bwt: Cannot open the ext SA file: %s!\n", extSAFname);
		fflush(stdout);
		exit(1);
	}
	bwtint_t sa0_index;
	bwtSA[0] = n;
	bwt_seq[0] = seq[bwtSA[0] - 1];
	// stream from file
	for (bwtint_t i = 1; i <= n; i++)
	{
		bwtint_t sa_i;
		if (fread(&sa_i, 5, 1, saFile) < 1)
			load_ext_sa_error(extSAFname); // 40-bit
		if (i % SA_INTERVAL == 0)
		{
			bwtSA[i / SA_INTERVAL] = sa_i;
		}
		// set SA to BWT
		if (sa_i == 0)
		{
			sa0_index = i;
			bwt_seq[i] = nt16_table['$'];
		}
		else
		{
			bwt_seq[i] = seq[sa_i - 1];
		}
	}
	return sa0_index;
}

// Construct the BWT index
bwt_t *construct_bwt(unsigned char *ref, const bwtint_t length, const char *extSAFname)
{
	bwt_t *BWT = (bwt_t *)calloc(1, sizeof(bwt_t));
	BWT->length = length + 1;

	// alloc suffix array
	BWT->num_sa = ceil(((double)BWT->length) / SA_INTERVAL);
	BWT->SA = (bwtint_t *)calloc(BWT->num_sa, sizeof(bwtint_t));
	if (BWT->SA == 0)
	{
		printf("construct_bwt: Could not allocate memory for the compressed SA, alloc'ing %" PRIbwtint_t " Mb\n", BWT->num_sa * sizeof(bwtint_t) / 1024 / 1024);
		fflush(stdout);
		exit(1);
	}
	// compute the BWT (and SA) of the sequence (output in ref)
	clock_t t = clock();
	if (extSAFname == NULL)
	{
		BWT->sa0_index = is_bwt(ref, length, BWT->SA);
		printf("SAIS time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	else
	{
		printf("Computing BWT from precomputed eSAIS SA \n");
		unsigned char *bwt_seq = (unsigned char *)malloc((length + 1) * sizeof(unsigned char));
		BWT->sa0_index = esa2bwt(ref, bwt_seq, length, BWT->SA, extSAFname);
		free(ref);
		ref = bwt_seq;
	}
	if (BWT->sa0_index < 0)
	{
		printf("SA construction failed\n");
		fflush(stdout);
		exit(1);
	}

	// compressed BWT code
	BWT->num_words = ceil(((double)BWT->length) / CHARS_PER_WORD);
	BWT->bwt = (uint32_t *)calloc(BWT->num_words, sizeof(uint32_t));

	// occurrence array
	BWT->num_occ = ceil(((double)BWT->length) / OCC_INTERVAL);
	BWT->O = (bwtint_t *)calloc(BWT->num_occ * ALPHABET_SIZE, sizeof(bwtint_t));

	if ((BWT->bwt == 0) || (BWT->O == 0) || (BWT->SA == 0))
	{
		printf("construct_bwt: Could not allocate memory for the BWT index, alloc'ing %" PRIbwtint_t " Mb\n", (BWT->num_words * sizeof(uint32_t) + BWT->num_occ * ALPHABET_SIZE * sizeof(bwtint_t)) / 1024 / 1024);
		fflush(stdout);
		exit(1);
	}

	// compress the BWT
	printf("Compressing BWT\n");
	t = clock();
	pack_word(ref, BWT->bwt, BWT->length);
	printf("BWT compression time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// compute the auxiliary C and O values
	printf("Computing C-array\n");
	t = clock();
	compute_C(BWT, ref);
	printf("C-array construction time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	free(ref);

	printf("Computing O-array\n");
	t = clock();
	compute_O(BWT);
	printf("O-array construction time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	return BWT;
}

void free_bwt(bwt_t *BWT)
{
	free(BWT->bwt);
	//_mm_free(BWT->bwt);
	free(BWT->O);
	if (BWT->SA)
		free(BWT->SA);
	free(BWT);
}

void print_bwt(const bwt_t *BWT)
{
	printf("BWT: \n");
	printf("length = %" PRIbwtint_t "\n", BWT->length);
	printf("num_words = %" PRIbwtint_t "\n", BWT->num_words);
	for (bwtint_t i = 0; i < BWT->length; i++)
	{
		if (i % 10000 != 0)
			continue;
		if (i == BWT->sa0_index)
		{
			printf("B[%" PRIbwtint_t "] = %c*\n", i, iupacChar[B(BWT, i)]);
		}
		else
		{
			printf("B[%" PRIbwtint_t "] = %c\n", i, iupacChar[B(BWT, i)]);
		}
	}
	printf("-----\n");
	for (bwtint_t i = 0; i < ALPHABET_SIZE; i++)
	{
		printf("C[%c] = %" PRIbwtint_t "\n", iupacChar[i], BWT->C[i]);
	}
	printf("-----\n");
	printf("num_occ = %" PRIbwtint_t "\n", BWT->num_occ);
	for (bwtint_t i = 0; i < BWT->num_occ; i++)
	{
		if (i % 10000 != 0)
			continue;
		printf("i = %" PRIbwtint_t ": ", i);
		for (bwtint_t j = 0; j < ALPHABET_SIZE; j++)
		{
			printf("O[%" PRIbwtint_t "] = %" PRIbwtint_t " ", j, BWT->O[i * ALPHABET_SIZE + j]);
		}
		printf("\n");
	}
	printf("-----\n");
	printf("num_sa = %" PRIbwtint_t "\n", BWT->num_sa);
	if (BWT->SA)
	{
		for (bwtint_t i = 0; i < BWT->num_sa; i++)
		{
			if (i % 10000 != 0)
				continue;
			printf("SA[%" PRIbwtint_t "] = %" PRIbwtint_t "\n", i, BWT->SA[i]);
		}
		printf("-----\n");
	}
}

// Count values
void compute_C(bwt_t *BWT, const unsigned char *seq)
{
	// count # occurrences of each character
	for (bwtint_t i = 0; i < BWT->length; i++)
	{
		if (i != BWT->sa0_index)
		{
			BWT->C[seq[i] + 1]++;
		}
	}
	// prefix sum
	for (bwtint_t i = 1; i <= ALPHABET_SIZE; i++)
	{
		BWT->C[i] += BWT->C[i - 1];
	}
}

// Occurrence values
void compute_O(bwt_t *BWT)
{
	bwtint_t occ[ALPHABET_SIZE] = {0};
	for (bwtint_t i = 0; i < BWT->length; i++)
	{
		const unsigned char c = B(BWT, i);
		if (i != BWT->sa0_index)
		{
			occ[c]++;
		}
		if (i % OCC_INTERVAL == 0)
		{
			memcpy(&(BWT->O[(i / OCC_INTERVAL) * ALPHABET_SIZE]), occ, ALPHABET_SIZE * sizeof(bwtint_t));
		}
	}
}

// Suffix Array values
void compute_SA(bwt_t *BWT)
{
	bwtint_t isa = 0;
	bwtint_t sa = BWT->length;
	for (bwtint_t i = 0; i < BWT->length; i++)
	{
		// apply invPsi until we get to a factor of SA_INTERVAL
		if (isa % SA_INTERVAL == 0)
		{
			BWT->SA[isa / SA_INTERVAL] = sa;
		}
		sa--;
		isa = invPsi(BWT, isa);
	}
	if (isa % SA_INTERVAL == 0)
	{
		BWT->SA[isa / SA_INTERVAL] = sa;
	}
}

// invPsi(i) = C(B[i]) + O(B[i],i)
bwtint_t invPsi(const bwt_t *BWT, const bwtint_t i)
{
	if (i == BWT->sa0_index)
	{
		return 0;
	}
	const unsigned char c = B(BWT, i);
	return C(BWT, c) + O(BWT, c, i);
}

// SA[i] = get the start position of the ith smallest suffix
bwtint_t SA(const bwt_t *BWT, bwtint_t i)
{
	// repeatedly apply invPsi until invPsi(i) is a factor of SA_INTERVAL
	bwtint_t j = 0;
	while (!(i % SA_INTERVAL == 0))
	{
		i = invPsi(BWT, i);
		j++;
	}
	// SA[i] = SA[invPsi(i) applied j times] + j
	return (BWT->SA[i / SA_INTERVAL] + j) % (BWT->length);
}

// C[c] := get the number of symbols lexicographically smaller than c
bwtint_t C(const bwt_t *BWT, const unsigned char c)
{
	return BWT->C[c];
}

// B[i] := get the char at index i in the compressed BWT
unsigned char B(const bwt_t *BWT, const bwtint_t i)
{
	const bwtint_t wordIndex = i / CHARS_PER_WORD;
	const bwtint_t charIndex = i - wordIndex * CHARS_PER_WORD;

	uint32_t w = BWT->bwt[wordIndex];
	w <<= charIndex * BITS_PER_CHAR;
	const unsigned char c = w >> (BITS_IN_WORD - BITS_PER_CHAR);
	return c;
}

// O(B[i], i) := get the occurrence count of char c at position i
inline bwtint_t O(const bwt_t *BWT, const unsigned char c, const bwtint_t i)
{
	if (i == BWT->length - 1)
	{
		return BWT->C[c + 1] - BWT->C[c];
	}
	if (i == (bwtint_t)-1)
	{
		return 0;
	}

	// closest index to i for which we have the O value
	const bwtint_t k = i / OCC_INTERVAL;
	bwtint_t o_c = BWT->O[k * ALPHABET_SIZE + c];

	if (c != 0)
	{ //&& ((k*OCC_INTERVAL) % CHARS_PER_WORD == 0)) {
		o_c += get_occ_count_opt(BWT, c, k * OCC_INTERVAL, i);
	}
	else
	{
		for (bwtint_t j = k * OCC_INTERVAL + 1; j <= i; j++)
		{
			if (j == BWT->sa0_index)
				continue;
			const unsigned char c_j = B(BWT, j);
			if (c_j == c)
			{
				o_c++;
			}
		}
	}
	return o_c;
}

void O_alphabet(const bwt_t *BWT, const bwtint_t i, int alphabet_size, bwtint_t *occ, int inc)
{
	if (i == BWT->length - 1)
	{
		occ[1] = BWT->C[2] + inc;
		occ[2] = BWT->C[3] + inc;
		occ[3] = BWT->C[4] + inc;
		occ[4] = BWT->C[5] + inc;
		occ[5] = BWT->C[6] + inc;
		occ[6] = BWT->C[7] + inc;
		occ[7] = BWT->C[8] + inc;
		occ[8] = BWT->C[9] + inc;
		occ[9] = BWT->C[10] + inc;
		occ[10] = BWT->C[11] + inc;
		occ[11] = BWT->C[12] + inc;
		occ[12] = BWT->C[13] + inc;
		occ[13] = BWT->C[14] + inc;
		occ[14] = BWT->C[15] + inc;
		occ[15] = BWT->C[16] + inc;
		return;
	}
	if (i == (bwtint_t)-1)
	{
		occ[1] = BWT->C[1] + inc;
		occ[2] = BWT->C[2] + inc;
		occ[3] = BWT->C[3] + inc;
		occ[4] = BWT->C[4] + inc;
		occ[5] = BWT->C[5] + inc;
		occ[6] = BWT->C[6] + inc;
		occ[7] = BWT->C[7] + inc;
		occ[8] = BWT->C[8] + inc;
		occ[9] = BWT->C[9] + inc;
		occ[10] = BWT->C[10] + inc;
		occ[11] = BWT->C[11] + inc;
		occ[12] = BWT->C[12] + inc;
		occ[13] = BWT->C[13] + inc;
		occ[14] = BWT->C[14] + inc;
		occ[15] = BWT->C[15] + inc;
		return;
	}

	// closest index to i for which we have the O value
	const bwtint_t k = i / OCC_INTERVAL;
	const bwtint_t s = k * ALPHABET_SIZE;
	/*for(int j = 1; j < alphabet_size; j++) {
		if((j == 5) || (j == 9) || (j == 13) || (j == 11)) {
			occ[j] = inc;
		} else {
			occ[j] = BWT->C[j] + BWT->O[s + j] + inc;
		}
	}*/
	get_occ_count_alphabet(BWT, k * OCC_INTERVAL, i, alphabet_size, occ);
	occ[1] += BWT->C[1] + BWT->O[s + 1] + inc;
	occ[2] += BWT->C[2] + BWT->O[s + 2] + inc;
	occ[3] += BWT->C[3] + BWT->O[s + 3] + inc;
	occ[4] += BWT->C[4] + BWT->O[s + 4] + inc;
	occ[5] += BWT->C[5] + inc;
	occ[6] += BWT->C[6] + BWT->O[s + 6] + inc;
	occ[7] += BWT->C[7] + BWT->O[s + 7] + inc;
	occ[8] += BWT->C[8] + BWT->O[s + 8] + inc;
	occ[9] += BWT->C[9] + inc;
	occ[10] += BWT->C[10] + BWT->O[s + 10] + inc;
	occ[11] += BWT->C[11] + inc;
	occ[12] += BWT->C[12] + BWT->O[s + 12] + inc;
	occ[13] += BWT->C[13] + inc;
	occ[14] += BWT->C[14] + BWT->O[s + 14] + inc;
	occ[15] += BWT->C[15] + BWT->O[s + 15] + inc;
}

void O_actg_alphabet(const bwt_t *BWT, const bwtint_t i, bwtint_t *occ, int inc)
{
	if (i == BWT->length - 1)
	{
		occ[1] = BWT->C[16] + inc;
		occ[2] = BWT->C[4] + inc;
		occ[3] = BWT->C[8] + inc;
		occ[4] = BWT->C[2] + inc;
		return;
	}
	if (i == (bwtint_t)-1)
	{
		occ[1] = BWT->C[15] + inc;
		occ[2] = BWT->C[3] + inc;
		occ[3] = BWT->C[7] + inc;
		occ[4] = BWT->C[1] + inc;
		return;
	}
	// closest index to i for which we have the O value
	const bwtint_t k = i / OCC_INTERVAL;
	const bwtint_t s = k * ALPHABET_SIZE;
	get_occ_count_actg_alphabet(BWT, k * OCC_INTERVAL, i, occ);
	occ[1] += BWT->C[15] + BWT->O[s + 15] + inc;
	occ[2] += BWT->C[3] + BWT->O[s + 3] + inc;
	occ[3] += BWT->C[7] + BWT->O[s + 7] + inc;
	occ[4] += BWT->C[1] + BWT->O[s + 1] + inc;
}

// L <= U
void O_LU(const bwt_t *BWT, const unsigned char c, const bwtint_t L, const bwtint_t U, bwtint_t *occL, bwtint_t *occU)
{

	if (L == (bwtint_t)-1)
	{
		*occL = 0;
		*occU = O(BWT, c, U);
		return;
	}
	if (U == BWT->length - 1)
	{
		*occL = O(BWT, c, L);
		*occU = BWT->C[c + 1] - BWT->C[c];
		return;
	}

	const bwtint_t k = L / OCC_INTERVAL;
	bwtint_t o_c = BWT->O[k * ALPHABET_SIZE + c]; // stored O value

	const bwtint_t start_index = k * OCC_INTERVAL;
	const bwtint_t bwt_start = start_index >> 3;
	const bwtint_t bwt_num_words = (L - start_index + 1) >> 3;

	// common L and U words
	for (bwtint_t i = bwt_start; i < bwt_start + bwt_num_words; i++)
	{
		const uint32_t w = BWT->bwt[i] ^ char_masks[c];
		o_c += BWT->occ_count_table[w & 0x0000FFFF];
		o_c += BWT->occ_count_table[w >> 16];
	}

	bwtint_t occ_l = o_c;
	bwtint_t occ_u = o_c;

	// partial word with L
	int num_chars = L - (start_index + (bwt_num_words << 3)) + 1;
	if (num_chars > 0)
	{
		// 0s where c
		const uint32_t w = (BWT->bwt[bwt_start + bwt_num_words] ^ char_masks[c]) | word_fractions_masks[num_chars - 1];
		occ_l += BWT->occ_count_table[w & 0x0000FFFF];
		occ_l += BWT->occ_count_table[w >> 16];
	}
	if ((BWT->bwt[bwt_start] >> 28) == c)
	{
		occ_l--;
	}
	*occL = occ_l;

	// remaining positions for U
	const bwtint_t next_word_pos = start_index + bwt_num_words * CHARS_PER_WORD;
	occ_u += get_occ_count_opt(BWT, c, next_word_pos, U);

	if (bwt_num_words != 0)
	{
		if ((BWT->bwt[bwt_start + bwt_num_words] >> 28) == c)
		{
			occ_u++;
		}
	}
	*occU = occ_u;
}

/* Optimizations */

// for all possible combinations of 2 bytes (total of 4 chars)
// compute the count of 0 in these 2 bytes
void generate_occ_table(bwt_t *BWT)
{
	char zero = 0;
	// for each possible combination of 2 bytes (4 chars)
	for (int i = 0; i < (1 << (2 * BITS_IN_BYTE)); i++)
	{
		uint8_t zero_count = 0;
		// count how many times zero occurs in i
		for (int j = 0; j < (2 * BITS_IN_BYTE) / BITS_PER_CHAR; j++)
		{
			zero_count += (zero == ((i >> (j * BITS_PER_CHAR)) & CHAR_MASK));
		}
		BWT->occ_count_table[i] = zero_count;
	}
}

// get the occurrence of character c in the interval (start_index, end_index]
// start_index is guaranteed to be the first position in the compressed BWT word
// and a multiple of OCC_INTERVAL
inline bwtint_t get_occ_count(const bwt_t *BWT, const unsigned char c, const bwtint_t start_index, const bwtint_t end_index)
{
	// the first word in the compressed BWT
	const bwtint_t bwt_start = start_index / CHARS_PER_WORD;
	// number of whole words of the compressed BWT that are in the interval
	const bwtint_t bwt_num_words = (end_index - start_index + 1) / CHARS_PER_WORD;

	bwtint_t total_count = 0;
	for (bwtint_t i = bwt_start; i < bwt_start + bwt_num_words; i++)
	{
		// xor the word with the mask for the character
		const uint32_t bwtword_xor_char = BWT->bwt[i] ^ char_masks[c];
		const uint32_t ls_key = bwtword_xor_char & LS_2BYTES_MASK;
		const uint32_t ms_key = bwtword_xor_char >> 16;
		total_count += BWT->occ_count_table[ls_key];
		total_count += BWT->occ_count_table[ms_key];
	}

	// count the remaining positions individually
	int num_chars = end_index - (start_index + (bwt_num_words << 3)) + 1;
	uint32_t w = BWT->bwt[bwt_start + bwt_num_words];
	for (int i = 0; i < num_chars; i++)
	{
		uint32_t w_i = w << (i * BITS_PER_CHAR);
		char c_i = w_i >> 28;
		if (c_i == c)
		{
			total_count++;
		}
	}

	if ((BWT->bwt[bwt_start] >> (BITS_IN_WORD - BITS_PER_CHAR)) == c)
	{
		total_count--;
	}

	return total_count;
}

bwtint_t get_occ_count_opt(const bwt_t *BWT, const unsigned char c, const bwtint_t start_index, const bwtint_t end_index)
{
	//_mm_prefetch((const char*)(BWT->bwt + start_index / CHARS_PER_WORD), _MM_HINT_NTA);

	const bwtint_t bwt_start = start_index >> 3;
	const bwtint_t bwt_num_words = (end_index - start_index + 1) >> 3;
	const int num_chars = end_index - (start_index + (bwt_num_words << 3)) + 1;

	bwtint_t total_count = 0;
	for (bwtint_t i = bwt_start; i < bwt_start + bwt_num_words; i++)
	{
		const uint32_t w = BWT->bwt[i] ^ char_masks[c];
		total_count += BWT->occ_count_table[w & 0x0000FFFF];
		total_count += BWT->occ_count_table[w >> 16];
	}

	if (num_chars > 0)
	{
		// 0s where c
		const uint32_t w = (BWT->bwt[bwt_start + bwt_num_words] ^ char_masks[c]) | word_fractions_masks[num_chars - 1];
		total_count += BWT->occ_count_table[w & 0x0000FFFF];
		total_count += BWT->occ_count_table[w >> 16];
	}

	if ((BWT->bwt[bwt_start] >> 28) == c)
	{
		total_count--;
	}
	return total_count;
}

bwtint_t get_occ_count_opt_sse(const bwt_t *BWT, const unsigned char c, const bwtint_t start_index, const bwtint_t end_index)
{
	// pre-fetch the compressed BWT interval (hint not to keep in lower-level caches)
	// loads 1 cache line (64 bytes? = 128 chars)
	_mm_prefetch((const char *)(BWT->bwt + start_index / CHARS_PER_WORD), _MM_HINT_NTA);

	const bwtint_t bwt_w0 = start_index / CHARS_PER_WORD;
	const bwtint_t num_words = (end_index - start_index + 1 + CHARS_PER_WORD - 1) / CHARS_PER_WORD;
	const bwtint_t num_128bits = (end_index - start_index + 1) / CHARS_PER_128BITS;
	const int num_chars = end_index - (start_index + (num_128bits * CHARS_PER_128BITS)) + 1;

	// SSE registers (128 * 4-bits/char = 16 32-bit words => 4 128-bit SSE registers)
	__m128i bwtR[4];
	__m128i xor_maskR = _mm_load_si128((__m128i *)(&char_masks128[c]));
	__m128i numchar_maskR = _mm_load_si128((__m128i *)(&word_fractions_masks128[num_chars - 1]));

	// XOR with the character mask results
	uint32_t w[16] __attribute__((aligned(16)));

	// full 128-bit batches
	for (int i = 0; i < num_128bits; i++)
	{
		bwtR[i] = _mm_load_si128((__m128i *)(BWT->bwt + bwt_w0 + i * 4));
		bwtR[i] = _mm_xor_si128(bwtR[i], xor_maskR);
		_mm_store_si128((__m128i *)&w[i * 4], bwtR[i]);
	}
	// last (partial) 128-bit batch
	if (num_chars > 0)
	{
		bwtR[num_128bits] = _mm_load_si128((__m128i *)(BWT->bwt + bwt_w0 + num_128bits * 4));
		bwtR[num_128bits] = _mm_xor_si128(bwtR[num_128bits], xor_maskR);
		bwtR[num_128bits] = _mm_or_si128(bwtR[num_128bits], numchar_maskR);
		_mm_store_si128((__m128i *)&w[num_128bits * 4], bwtR[num_128bits]);
	}

	bwtint_t total_count = 0;
	for (int j = 0; j < num_words; j++)
	{
		total_count += BWT->occ_count_table[w[j] & 0x0000FFFF];
		total_count += BWT->occ_count_table[w[j] >> 16];
	}

	// bwt_w0[0] is counted twice
	if ((BWT->bwt[bwt_w0] >> 28) == c)
	{
		total_count--;
	}
	return total_count;
}

void get_occ_count_actg_alphabet(const bwt_t *BWT, const bwtint_t start_index, const bwtint_t end_index, bwtint_t *occ)
{
	// the first word in the compressed BWT
	const bwtint_t bwt_start = start_index >> 3;
	// number of whole words of the compressed BWT that are in the interval
	const bwtint_t bwt_num_words = (end_index - start_index + 1) >> 3;

	occ[(BWT->bwt[bwt_start] >> 28)]--; // position already counted in O array
	// map to A/C/T/G alphabet
	occ[4] = occ[1];
	occ[1] = occ[15];
	occ[2] = occ[3];
	occ[3] = occ[7];

	// full compressed words
	for (bwtint_t i = bwt_start; i < bwt_start + bwt_num_words; i++)
	{
		const uint32_t w = BWT->bwt[i];
		const uint32_t wA = w ^ CHAR_MASK_A;
		const uint32_t wG = w ^ CHAR_MASK_G;
		const uint32_t wC = w ^ CHAR_MASK_C;
		const uint32_t wT = w ^ CHAR_MASK_T;
		occ[1] += BWT->occ_count_table[wA & 0x0000FFFF] + BWT->occ_count_table[wA >> 16];
		occ[2] += BWT->occ_count_table[wG & 0x0000FFFF] + BWT->occ_count_table[wG >> 16];
		occ[3] += BWT->occ_count_table[wC & 0x0000FFFF] + BWT->occ_count_table[wC >> 16];
		occ[4] += BWT->occ_count_table[wT & 0x0000FFFF] + BWT->occ_count_table[wT >> 16];
	}

	// last word (could be partially full)
	int num_chars = end_index - (start_index + (bwt_num_words << 3)) + 1;
	if (num_chars > 0)
	{
		const uint32_t w = BWT->bwt[bwt_start + bwt_num_words];
		const uint32_t w_mask = word_fractions_masks[num_chars - 1];
		const uint32_t wA = (w ^ CHAR_MASK_A) | w_mask;
		const uint32_t wG = (w ^ CHAR_MASK_G) | w_mask;
		const uint32_t wC = (w ^ CHAR_MASK_C) | w_mask;
		const uint32_t wT = (w ^ CHAR_MASK_T) | w_mask;
		occ[1] += BWT->occ_count_table[wA & 0x0000FFFF] + BWT->occ_count_table[wA >> 16];
		occ[2] += BWT->occ_count_table[wG & 0x0000FFFF] + BWT->occ_count_table[wG >> 16];
		occ[3] += BWT->occ_count_table[wC & 0x0000FFFF] + BWT->occ_count_table[wC >> 16];
		occ[4] += BWT->occ_count_table[wT & 0x0000FFFF] + BWT->occ_count_table[wT >> 16];
	}
}

void get_occ_count_alphabet(const bwt_t *BWT, const bwtint_t start_index, const bwtint_t end_index, int alphabet_size, bwtint_t *occ)
{
	// the first word in the compressed BWT
	const bwtint_t bwt_start = start_index >> 3;
	// number of whole words of the compressed BWT that are in the interval
	const bwtint_t bwt_num_words = (end_index - start_index + 1) >> 3;

	for (bwtint_t i = bwt_start; i < bwt_start + bwt_num_words; i++)
	{
		// xor the word with the mask for the character
		const uint32_t w = BWT->bwt[i];
		const uint32_t w1 = w ^ char_masks[1];
		const uint32_t w2 = w ^ char_masks[2];
		const uint32_t w3 = w ^ char_masks[3];
		const uint32_t w4 = w ^ char_masks[4];
		//const uint32_t w5 = w ^ char_masks[5];
		const uint32_t w6 = w ^ char_masks[6];
		const uint32_t w7 = w ^ char_masks[7];
		const uint32_t w8 = w ^ char_masks[8];
		//const uint32_t w9 = w ^ char_masks[9];
		const uint32_t w10 = w ^ char_masks[10];
		//const uint32_t w11 = w ^ char_masks[11];
		const uint32_t w12 = w ^ char_masks[12];
		//const uint32_t w13 = w ^ char_masks[13];
		const uint32_t w14 = w ^ char_masks[14];
		const uint32_t w15 = w ^ char_masks[15];

		occ[1] += BWT->occ_count_table[w1 & 0x0000FFFF] + BWT->occ_count_table[w1 >> 16];
		occ[2] += BWT->occ_count_table[w2 & 0x0000FFFF] + BWT->occ_count_table[w2 >> 16];
		occ[3] += BWT->occ_count_table[w3 & 0x0000FFFF] + BWT->occ_count_table[w3 >> 16];
		occ[4] += BWT->occ_count_table[w4 & 0x0000FFFF] + BWT->occ_count_table[w4 >> 16];
		//occ[5] += BWT->occ_count_table[w5 & 0x0000FFFF] + BWT->occ_count_table[w5 >> 16];
		occ[6] += BWT->occ_count_table[w6 & 0x0000FFFF] + BWT->occ_count_table[w6 >> 16];
		occ[7] += BWT->occ_count_table[w7 & 0x0000FFFF] + BWT->occ_count_table[w7 >> 16];
		occ[8] += BWT->occ_count_table[w8 & 0x0000FFFF] + BWT->occ_count_table[w8 >> 16];
		//occ[9] += BWT->occ_count_table[w9 & 0x0000FFFF] + BWT->occ_count_table[w9 >> 16];
		occ[10] += BWT->occ_count_table[w10 & 0x0000FFFF] + BWT->occ_count_table[w10 >> 16];
		//occ[11] += BWT->occ_count_table[w11 & 0x0000FFFF] + BWT->occ_count_table[w11 >> 16];
		occ[12] += BWT->occ_count_table[w12 & 0x0000FFFF] + BWT->occ_count_table[w12 >> 16];
		//occ[13] += BWT->occ_count_table[w13 & 0x0000FFFF] + BWT->occ_count_table[w13 >> 16];
		occ[14] += BWT->occ_count_table[w14 & 0x0000FFFF] + BWT->occ_count_table[w14 >> 16];
		occ[15] += BWT->occ_count_table[w15 & 0x0000FFFF] + BWT->occ_count_table[w15 >> 16];

		/*for(int j = 1; j < alphabet_size; j++) {
			if((j == 5) || (j == 9) || (j == 13) || (j == 11)) continue;
			const uint32_t wchar = BWT->bwt[i] ^ char_masks[j];
			occ[j] += BWT->occ_count_table[wchar & 0x0000FFFF] + BWT->occ_count_table[wchar >> 16];
		}*/
	}

	int num_chars = end_index - (start_index + (bwt_num_words << 3)) + 1;
	if (num_chars > 0)
	{
		const uint32_t w = BWT->bwt[bwt_start + bwt_num_words];
		const uint32_t w_mask = word_fractions_masks[num_chars - 1];

		const uint32_t w1 = (w ^ char_masks[1]) | w_mask;
		const uint32_t w2 = (w ^ char_masks[2]) | w_mask;
		const uint32_t w3 = (w ^ char_masks[3]) | w_mask;
		const uint32_t w4 = (w ^ char_masks[4]) | w_mask;
		//const uint32_t w5 = (w ^ char_masks[5]) | w_mask;
		const uint32_t w6 = (w ^ char_masks[6]) | w_mask;
		const uint32_t w7 = (w ^ char_masks[7]) | w_mask;
		const uint32_t w8 = (w ^ char_masks[8]) | w_mask;
		//const uint32_t w9 = (w ^ char_masks[9]) | w_mask;
		const uint32_t w10 = (w ^ char_masks[10]) | w_mask;
		//const uint32_t w11 = (w ^ char_masks[11]) | w_mask;
		const uint32_t w12 = (w ^ char_masks[12]) | w_mask;
		//const uint32_t w13 = (w ^ char_masks[13]) | w_mask;
		const uint32_t w14 = (w ^ char_masks[14]) | w_mask;
		const uint32_t w15 = (w ^ char_masks[15]) | w_mask;

		occ[1] += BWT->occ_count_table[w1 & 0x0000FFFF] + BWT->occ_count_table[w1 >> 16];
		occ[2] += BWT->occ_count_table[w2 & 0x0000FFFF] + BWT->occ_count_table[w2 >> 16];
		occ[3] += BWT->occ_count_table[w3 & 0x0000FFFF] + BWT->occ_count_table[w3 >> 16];
		occ[4] += BWT->occ_count_table[w4 & 0x0000FFFF] + BWT->occ_count_table[w4 >> 16];
		//occ[5] += BWT->occ_count_table[w5 & 0x0000FFFF] + BWT->occ_count_table[w5 >> 16];
		occ[6] += BWT->occ_count_table[w6 & 0x0000FFFF] + BWT->occ_count_table[w6 >> 16];
		occ[7] += BWT->occ_count_table[w7 & 0x0000FFFF] + BWT->occ_count_table[w7 >> 16];
		occ[8] += BWT->occ_count_table[w8 & 0x0000FFFF] + BWT->occ_count_table[w8 >> 16];
		//occ[9] += BWT->occ_count_table[w9 & 0x0000FFFF] + BWT->occ_count_table[w9 >> 16];
		occ[10] += BWT->occ_count_table[w10 & 0x0000FFFF] + BWT->occ_count_table[w10 >> 16];
		//occ[11] += BWT->occ_count_table[w11 & 0x0000FFFF] + BWT->occ_count_table[w11 >> 16];
		occ[12] += BWT->occ_count_table[w12 & 0x0000FFFF] + BWT->occ_count_table[w12 >> 16];
		//occ[13] += BWT->occ_count_table[w13 & 0x0000FFFF] + BWT->occ_count_table[w13 >> 16];
		occ[14] += BWT->occ_count_table[w14 & 0x0000FFFF] + BWT->occ_count_table[w14 >> 16];
		occ[15] += BWT->occ_count_table[w15 & 0x0000FFFF] + BWT->occ_count_table[w15 >> 16];
		/*for(int j = 1; j < alphabet_size; j++) {
			if((j == 5) || (j == 9) || (j == 13) || (j == 11)) continue;
			const uint32_t wchar = (w ^ char_masks[j]) | w_mask;
			occ[j] += BWT->occ_count_table[wchar & 0x0000FFFF] + BWT->occ_count_table[wchar >> 16];
		}*/
	}

	occ[(BWT->bwt[bwt_start] >> 28)]--;
}
