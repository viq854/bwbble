/*
 * sais.c for sais-lite
 * Copyright (c) 2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdlib.h>
#include "bwt.h"
#include "io.h"

typedef unsigned char ubyte_t;
typedef long int seqint_t;

#define chr(i) (cs == sizeof(seqint_t) ? ((const seqint_t *)T)[i]:((const unsigned char *)T)[i])

/* find the start or end of each bucket */
static void getCounts(const unsigned char *T, seqint_t *C, seqint_t n, seqint_t k, int cs)
{
	seqint_t i;
	for (i = 0; i < k; ++i) C[i] = 0;
	for (i = 0; i < n; ++i) ++C[chr(i)];
}
static void getBuckets(const seqint_t *C, seqint_t *B, seqint_t k, seqint_t end)
{
	seqint_t i, sum = 0;
	if (end) {
		for (i = 0; i < k; ++i) {
			sum += C[i];
			B[i] = sum;
		}
	} else {
		for (i = 0; i < k; ++i) {
			sum += C[i];
			B[i] = sum - C[i];
		}
	}
}

/* compute SA */
static void induceSA(const unsigned char *T, seqint_t *SA, seqint_t *C, seqint_t *B, seqint_t n, seqint_t k, int cs)
{
	seqint_t *b;//, i, j;
	seqint_t  c0, c1;
	seqint_t i;
	seqint_t j;
	/* compute SAl */
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 0);	/* find starts of buckets */
	j = n - 1;
	
        b = SA + B[c1 = chr(j)];
	*b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
	for (i = 0; i < n; ++i) {
		j = SA[i], SA[i] = ~j;
		if (0 < j) {
			--j;
			if ((c0 = chr(j)) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}
			*b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
		}
	}
	/* compute SAs */
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	/* find ends of buckets */
	
	//long int iL;
	for (i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
		if (0 < (j = SA[i])) {
			--j;
			if ((c0 = chr(j)) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}
			*--b = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
		} else SA[i] = ~j;
	}
}

/*
 * find the suffix array SA of T[0..n-1] in {0..k-1}^n use a working
 * space (excluding T and SA) of at most 2n+O(1) for a constant alphabet
 */
static int sais_main(const unsigned char *T, seqint_t *SA, seqint_t fs, seqint_t n, seqint_t k, int cs)
{
	seqint_t *C, *B, *RA;
	seqint_t  j, c, m, p, q, plen, qlen, name;
	long int i;
	seqint_t  c0, c1;
	seqint_t  diff;

	/* stage 1: reduce the problem by at least 1/2 sort all the
	 * S-substrings */
	if (k <= fs) {
		C = SA + n;
		B = (k <= (fs - k)) ? C + k : C;
	} else if ((C = B = (seqint_t *) malloc(k * sizeof(seqint_t))) == NULL) return -2;
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	/* find ends of buckets */
	for (i = 0; i < n; ++i) SA[i] = 0;
	for (i = n - 2, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr(i)) < (c1 + c)) c = 1;
		else if (c != 0) SA[--B[c1]] = i + 1, c = 0;
	}
	induceSA(T, SA, C, B, n, k, cs);
	if (fs < k) free(C);
	/* compact all the sorted substrings into the first m items of SA
	 * 2*m must be not larger than n (proveable) */
	for (i = 0, m = 0; i < n; ++i) {
		p = SA[i];
		if ((0 < p) && (chr(p - 1) > (c0 = chr(p)))) {
			for (j = p + 1; (j < n) && (c0 == (c1 = chr(j))); ++j);
			if ((j < n) && (c0 < c1)) SA[m++] = p;
		}
	}
	for (i = m; i < n; ++i) SA[i] = 0;	/* init the name array buffer */
	/* store the length of all substrings */
	for (i = n - 2, j = n, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr(i)) < (c1 + c)) c = 1;
		else if (c != 0) {
			SA[m + ((i + 1) >> 1)] = j - i - 1;
			j = i + 1;
			c = 0;
		}
	}
	/* find the lexicographic names of all substrings */
	for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
		p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
		if (plen == qlen) {
			for (j = 0; (j < plen) && (chr(p + j) == chr(q + j)); j++);
			if (j == plen) diff = 0;
		}
		if (diff != 0) ++name, q = p, qlen = plen;
		SA[m + (p >> 1)] = name;
	}

	/* stage 2: solve the reduced problem recurse if names are not yet
	 * unique */
	if (name < m) {
		RA = SA + n + fs - m;
		for (i = n - 1, j = m - 1; m <= i; --i) {
			if (SA[i] != 0) RA[j--] = SA[i] - 1;
		}
		if (sais_main((unsigned char *) RA, SA, fs + n - m * 2, m, name, sizeof(seqint_t)) != 0) return -2;
		for (i = n - 2, j = m - 1, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
			if ((c0 = chr(i)) < (c1 + c)) c = 1;
			else if (c != 0) RA[j--] = i + 1, c = 0; /* get p1 */
		}
		for (i = 0; i < m; ++i) SA[i] = RA[SA[i]]; /* get index */
	}
	/* stage 3: induce the result for the original problem */
	if (k <= fs) {
		C = SA + n;
		B = (k <= (fs - k)) ? C + k : C;
	} else if ((C = B = (seqint_t *) malloc(k * sizeof(seqint_t))) == NULL) return -2;
	/* put all left-most S characters into their buckets */
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	/* find ends of buckets */
	for (i = m; i < n; ++i) SA[i] = 0; /* init SA[m..n-1] */
	for (i = m - 1; 0 <= i; --i) {
		j = SA[i], SA[i] = 0;
		SA[--B[chr(j)]] = j;
	}
	induceSA(T, SA, C, B, n, k, cs);
	if (fs < k) free(C);
	return 0;
}

/**
 * Constructs the suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n] The output array of suffixes.
 * @param n The length of the given string.
 * @return 0 if no error occurred
 */
int is_sa(const ubyte_t *T, seqint_t *SA, seqint_t n)
{
	if ((T == NULL) || (SA == NULL) || (n < 0)) return -1;
	SA[0] = n;
	if (n <= 1) {
		if (n == 1) SA[1] = 0;
		return 0;
	}
	return sais_main(T, SA+1, 0, n, 256, 1);
}

/**
 * Constructs the burrows-wheeler transformed string of a given string.
 * @param T[0..n-1] The input string.
 * @param n The length of the given string.
 * @return The primary index if no error occurred, -1 or -2 otherwise.
 */
bwtint_t is_bwt(unsigned char *T, bwtint_t n, bwtint_t* bwtSA) {
	seqint_t *SA = (seqint_t*) calloc(n+1, sizeof(seqint_t));
	if(SA == 0) {
		printf("is_bwt: Could not allocate memory for the BWT index construction (in SAIS), alloc'ing %" PRIbwtint_t " Mb\n", (n+1)*sizeof(seqint_t)/1024/1024);
		exit(1);
	}
	is_sa(T, SA, n);

	bwtint_t sa0_index;
	for (bwtint_t i = 0; i <= n; i++) {
		if(i % SA_INTERVAL == 0) {
			bwtSA[i/SA_INTERVAL] = SA[i];
		}
		// set SA to BWT
		if (SA[i] == 0) {
			sa0_index = i;
			SA[i] = nt16_table['$'];
		}
		else {
			SA[i] = T[SA[i] - 1];
		}
	}

	// set T to BWT
	for (bwtint_t i = 0; i <= n; i++) {
		T[i] = SA[i];
	}
	free(SA);
	return sa0_index;
}

