//-------------------------------------------------------
// BWT Indexing
// Victoria Popic (viq@stanford.edu), Apr 2012
//-------------------------------------------------------

#ifndef BWTSNP_BWT_H
#define BWTSNP_BWT_H

#include <stdint.h>
#include "common.h"
#include "io.h"

// Occurrence values interval (only values O(_,k) where k is a factor of OCC_INTERVAL are stored)
#define OCC_INTERVAL 128
// Suffix array values interval (SA values will be stored only at this interval)
#define SA_INTERVAL 32

// BWT
typedef struct {
	// BWT length = reference length + 1 (in chars)
	bwtint_t length;
	// compressed BWT length (in words)
	bwtint_t num_words;
	// compressed BWT (each BITS_PER_CHAR bits correspond to a character in the sequence)
	uint32_t *bwt;
	// C(a) := # symbols smaller than a
	bwtint_t C[ALPHABET_SIZE + 1];
	// O(a,i) := # occurrences of a in BWT[0,i]
	bwtint_t* O;
	// Number of O values stored
	bwtint_t num_occ;
	// Character occurrence count table
	uint8_t occ_count_table[(1 << (BITS_IN_BYTE << 1))];
	// SA values
	bwtint_t *SA;
	// Number of SA values stored
	bwtint_t num_sa;
	// Index i s.t. SA[i] = 0
	bwtint_t sa0_index;
} bwt_t;

int index_bwt(const char *fastaFname, const char* extSAFname);
bwt_t* construct_bwt(unsigned char *seq, const bwtint_t length, const char* extSAFname);
void free_bwt(bwt_t* BWT);
void store_bwt(const bwt_t* BWT, const char* bwtFname);
bwt_t* load_bwt(const char* bwtFname, const int loadSA);
void print_bwt(const bwt_t *BWT);
bwtint_t is_bwt(unsigned char *T, bwtint_t n, bwtint_t* SA); // IS library method for constructing the transform

unsigned char B(const bwt_t* BWT, const bwtint_t i);
bwtint_t O(const bwt_t* BWT, const unsigned char c, bwtint_t i);
void O_alphabet(const bwt_t* BWT, const bwtint_t i, int alphabet_size, bwtint_t* occ, int L);
void O_actg_alphabet(const bwt_t* BWT, const bwtint_t i, bwtint_t* occ, int L);
void O_LU(const bwt_t* BWT, const unsigned char c, const bwtint_t L, const bwtint_t U, bwtint_t* occL, bwtint_t* occU);
bwtint_t C(const bwt_t* BWT, const unsigned char c);
bwtint_t SA(const bwt_t* BWT, const bwtint_t i);
bwtint_t invPsi(const bwt_t* BWT, const bwtint_t i);


#endif
