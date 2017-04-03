/*
HattCI is a statistical program used to identify attC sites of DNA sequences.
Copyright 2015 Mikael Wallroth, Mariana Buongermino Pereira

This file is part of HattCI.

HattCI is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HattCI is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HattCI.  If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------

Wrapper for processing sequences with viterbi and calcResults.
Returns when all sequences in seqs has been processed.

Parameters:
seqs 		(input)         char array of length seqLen
S			(input)         int, number of states
est_ini		(input)         double array of length S, initial probabilities
est_emi		(input)         double array of arrays, emission probabilities
est_trans	(input)         double array, transition probabilities
f			(input)         double array of arrays,
Dmin		(input)         double array, min length of each state
Dmax		(input)         double array, max length of each state
T		    (input)         int, length of each sequence
nseqs       (input)         int, number of sequences
headers     (input)         char vector of vectors, headers of each sequence
prevRead    (input/output)  int, number of seqs previously read
prevHits    (input/output)  int, number of hits previously found (original seq)
prevHitsRev (input/output)  int, number of hits previously found (reversed seq)
reverse     (input)         int, 1 if to check both strands, 0 else
output      (input)         char, output file to write results to
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void viterbi(int *, double *, char *, int, double *, double **, double *,
                double **, int *, int *, int);
int calcResults(int **, double **, char **, char **, int, int *, int, char **,
    int, int *, int *, int, char *, double **, double *, double *);

int hmm(char **seqs, int S, double *est_ini, double **est_emi, double *est_trans,
    double **f, int *Dmin, int *Dmax, int *T, int nseq, char **headers,
    int prevRead, int *prevHits, int *prevHitsRev, int reverse, char output[]) {

    int **resPaths, status, i;     // collection of paths for each sequence
    double **resLogV;   // collection of logV for each sequence
    char **revSeqs;

    // calculate the minimum required length to contain attC site
    int minLength = 0;
    for (int s = 0; s < S; s++) {
        minLength += Dmin[s];
    }

    if (reverse == 0) { // complementary strand is not to be run
        // allocate memory to store paths and logV for each sequence
        resPaths = malloc(nseq * sizeof(int *));
        resLogV = malloc(nseq * sizeof(double *));

        #pragma omp parallel for private(i) shared(nseq, T, seqs, resPaths, resLogV, S, est_ini, est_emi, est_trans, f, Dmin, Dmax)
        for (i = 0; i < nseq; i++) {
            resPaths[i] = malloc((T[i] + 1) * sizeof(int));
            resLogV[i] = malloc((T[i] + 1) * S * sizeof(double));

            // check if sequence is large enough to contain attC site
            if (T[i] >= minLength - 1) {
                // run viterbi for each sequence
                viterbi(resPaths[i], resLogV[i], seqs[i], S, est_ini, est_emi,
                    est_trans, f, Dmin, Dmax, T[i]);
            } else { // sequence less then minimum length, can't contain attC site
                for (int t = 0; t < T[i]; t++) {
                    resPaths[i][t] = 7;
                }
            }
        }
    } else { // complementary strand is to be run
        // allocate a vector to store the reversed sequence
        revSeqs = malloc (nseq * sizeof(char *));

        // allocate memory to store paths and logV for each sequence and its
        // complementary sequence (top/bottom strand)
        resPaths = malloc(2 * nseq * sizeof(int *));
        resLogV = malloc(2 * nseq * sizeof(double *));

        #pragma omp parallel for private(i) shared(nseq, T, seqs, resPaths, resLogV, S, est_ini, est_emi, est_trans, f, Dmin, Dmax)
        for (i = 0; i < nseq; i++) {
            resPaths[2 * i] = malloc((T[i] + 1) * sizeof(int));
            resLogV[2 * i] = malloc((T[i] + 1) * S * sizeof(double));
            resPaths[2 * i + 1] = malloc((T[i] + 1) * sizeof(int));
            resLogV[2 * i + 1] = malloc((T[i] + 1) * S * sizeof(double));
            revSeqs[i] = malloc(T[i] * sizeof(char));

            // check if sequence large enough to contain attC site
            if (T[i] >= minLength - 1) {
                // run viterbi for the sequence
                viterbi(resPaths[2 * i], resLogV[2 * i], seqs[i], S, est_ini,
                    est_emi, est_trans, f, Dmin, Dmax, T[i]);

                // transform the sequence to reversed complementary sequence
                for (int t = 0; t < T[i]; t++) {
                    if (seqs[i][T[i] - t - 1] == 'A' ||
                                            seqs[i][T[i] - t - 1] == 'a') {
                        revSeqs[i][t] = 'T';
                    } else if (seqs[i][T[i] - t - 1] == 'C' ||
                                            seqs[i][T[i] - t - 1] == 'c') {
                        revSeqs[i][t] = 'G';
                    } else if (seqs[i][T[i] - t - 1] == 'G' ||
                                            seqs[i][T[i] - t - 1] == 'g') {
                        revSeqs[i][t] = 'C';
                    } else {
                        revSeqs[i][t] = 'A';
                    }
                }

                // run viterbi for the reversed complementary sequence
                viterbi(resPaths[2 * i + 1], resLogV[2 * i + 1], revSeqs[i], S,
                    est_ini, est_emi, est_trans, f, Dmin, Dmax, T[i]);
            } else { // sequence less then minimum length, can't contain attC site
                for (int t = 0; t < T[i]; t++) {
                    resPaths[2 * i][t] = 7;
                    resPaths[2 * i + 1][t] = 7;
                }
            }
        }
    }

    // process the results
    if ((status = calcResults(resPaths, resLogV, seqs, revSeqs, nseq, T, S,
          headers, prevRead, prevHits, prevHitsRev, reverse, output, est_emi,
          est_trans, est_ini)) != 0) {
          printf("Error in calcResults!\n");
          return -1;
    }


    // free memory allocated in this function
    if (reverse == 1) {
        for (i = 0; i < nseq; i++) {
            free(revSeqs[i]);
        }
        free(revSeqs);
    }

    for (i = 0; i < (1 + reverse) * nseq; i++) {
        free(resPaths[i]);
        free(resLogV[i]);
    }
    free(resPaths);
    free(resLogV);

    return 0;

}
