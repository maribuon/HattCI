/*
HattCI is a statistical program used to identify attC sites of DNA sequences.
Copyright 2015 Mikael Wallroth

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

Process results from viterbi.
Prints properties of each attC site, such as difference in logV score and
coordinates of the states of the attC site.
Prints the sequence segments corresponding to each found attC site
together with properties of the attC sites.

Parameters:
resPaths    (input)         int array of arrays, resulting paths from viterbi
resLogV     (input)         double array of arrays, resulting logV from viterbi
seqs        (input)         char array of arrays, sequences to process
revSeqs     (input)         char array of arrays, reversed sequences to process
nseq        (input)         int, number of sequences
T           (input)         int array, length of each sequence
S           (input)         int, number of states
headers     (input)         char vector of vectors, headers of each sequence
prevRead    (input/output)  int, number of seqs previously read
totHits     (input/output)  int, input: # of hits prev found, output: # of hits
                                found including this iteration (original seq)
totHitsRev  (input/output)  int, input: # of hits prev found, output: # of hits
                                found including this iteration (reversed seq)
reverse     (input)         int, 1 if to check both strands, 0 else
output      (input)         char array, name of output file
est_emi     (input)         double array of arrays, emission log probabilities
est_trans	(input)         double array, transition log probabilities
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int calcResults(int **resPaths, double **resLogV, char **seqs, char **revSeqs,
    int nseq, int *T, int S, char **headers, int prevRead, int *totHits,
    int *totHitsRev, int reverse, char output[], double **est_emi,
    double *est_trans, double *est_ini) {


    int hits = 0, emi_rows[8] = {1, 4, 4, 16, 4, 4, 1, 4}, revHits = 0, *outN,
        *outIni, *outEnd, tt, *revIndex, *tmpVec;
    FILE *outfile;
    FILE *outTmp;
	 	FILE *outTmpFasta;

    // allocate memory to store results from found attC sites
    // the memory allocated is increased later if necessary
    outN = malloc(nseq * sizeof(int)); // which sequence
    outIni = malloc(nseq * sizeof(int)); // first element of the attC site
    outEnd = malloc(nseq * sizeof(int)); // last element of the attC sequence
    revIndex = malloc(nseq * sizeof(int)); // top or bottom strand (0 / 1)

    int size = nseq;
    int k = 0, t, i;

    // check for attC sites found in each sequence
    if (reverse == 0) { // no reversed
        for (i = 0; i < nseq; i++) {
            t = 1;
            // check if first element is R'' in an attC site
            if (resPaths[i][0] == 0) {
                hits++;

                // check if number of hits > currently allocated memory
                if (k >= size) {
                    size = size * 2;
                    if (!(tmpVec = realloc(outN, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outN = tmpVec;
                    if (!(tmpVec = realloc(outIni, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outIni = tmpVec;
                    if (!(tmpVec = realloc(outEnd, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outEnd = tmpVec;
                    if (!(tmpVec = realloc(revIndex, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    revIndex = tmpVec;
                }
                // store the sequence number and the first element of attC site
                // reversed is not true
                outN[k] = i + 1;
                outIni[k] = 1;
                revIndex[k] = 0;

                // find the last element of the attC site and store it
                for (tt = 0; tt < T[i] - 1; tt++) {
                    // the element that is in R' with the following in null state
                    if (resPaths[i][tt] == 6 && resPaths[i][tt+1] == 7) {
                        outEnd[k] = tt + 1;
                        break;
                    }
                    // the last element in the sequence is part of the attC site
                    if (tt == (T[i] - 2) && resPaths[i][T[i]-1] != 7) {
                        outEnd[k] = T[i];
                    }
                }
                k++;

                // set the current element to the first elemenet not in the attC site
                t = tt + 1;

            }
            // continue searching the sequence for attC sites, by looking for
            // R'' in the remaining elements
            while (t < T[i]) {
                // check if previous element is the null state and current element
                // is R'', if so, new attC site found
                if (resPaths[i][t-1] == 7 && resPaths[i][t] == 0) {
                    // increase memory if necessary
                    if (k >= size) {
                        size = size * 2;
                        if (!(tmpVec = realloc(outN, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outN = tmpVec;
                        if (!(tmpVec = realloc(outIni, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outIni = tmpVec;
                        if (!(tmpVec = realloc(outEnd, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outEnd = tmpVec;
                        if (!(tmpVec = realloc(revIndex, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        revIndex = tmpVec;
                    }
                    hits++;

                    // store sequence number and first element
                    outN[k] = i + 1;
                    outIni[k] = t + 1;
                    revIndex[k] = 0; // reversed = false

                    // find the last element of the attC site and store it
                    for (tt = t; tt < T[i] - 1; tt++) {
                        // the element that is in R' with the following in null state
                        if (resPaths[i][tt] == 6 && resPaths[i][tt+1] == 7) {
                            outEnd[k] = tt + 1;
                            break;
                        }
                        // the last element in the sequence is part of the attC site
                        if (tt == (T[i] - 2) && resPaths[i][T[i]-1] != 7) {
                            outEnd[k] = T[i];
                        }
                    }
                    k++;
                    // continue the search for new attC sites on the remaining
                    // elements in the sequence
                    t = tt + 1;
                } else { // element is not in R''
                    t++;
                }
            }
        }
    } else { // reversed, check two paths per iteration (original + complementary)
        int seqIndex;
        for (i = 0; i < nseq; i++) {
            seqIndex = 2 * i; // check the original sequence (even)
            t = 1;
            // check if first element is R'' in an attC site
            if (resPaths[seqIndex][0] == 0) {
                hits++;
                // check if number of hits > currently allocated memory
                if (k >= size) {
                    size = size * 2;
                    if (!(tmpVec = realloc(outN, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outN = tmpVec;
                    if (!(tmpVec = realloc(outIni, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outIni = tmpVec;
                    if (!(tmpVec = realloc(outEnd, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outEnd = tmpVec;
                    if (!(tmpVec = realloc(revIndex, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    revIndex = tmpVec;
                }
                // store the sequence number and the first element of attC site
                outN[k] = i + 1;
                outIni[k] = 1;
                revIndex[k] = 0; // reversed is not true

                // find the last element of the attC site and store it
                for (tt = 0; tt < T[i] - 1; tt++) {

                    // the element that is in R' with the following in null state
                    if (resPaths[seqIndex][tt] == 6 && resPaths[seqIndex][tt+1] == 7) {
                        outEnd[k] = tt + 1;
                        break;
                    }
                    // the last element in the sequence is part of the attC site
                    if (tt == (T[i] - 2) && resPaths[seqIndex][T[i]-1] != 7) {
                        outEnd[k] = T[i];
                    }
                }

                k++;
                // continue the search for new attC sites on the remaining
                // elements in the sequence
                t = tt + 1;
            }
            // continue searching the sequence for attC sites, by looking for
            // R'' in the remaining elements
            while (t < T[i]) {
                // check if previous element is the null state and current element
                // is R'', if so, new attC site found
                if (resPaths[seqIndex][t-1] == 7 && resPaths[seqIndex][t] == 0) {
                    // increase memory if necessary
                    if (k >= size) {
                        size = size * 2;
                        if (!(tmpVec = realloc(outN, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outN = tmpVec;
                        if (!(tmpVec = realloc(outIni, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outIni = tmpVec;
                        if (!(tmpVec = realloc(outEnd, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outEnd = tmpVec;
                        if (!(tmpVec = realloc(revIndex, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        revIndex = tmpVec;
                    }
                    hits++;
                    // store sequence number and first element
                    outN[k] = i + 1;
                    outIni[k] = t + 1;
                    revIndex[k] = 0; // reversed = false

                    // find the last element of the attC site and store it
                    for (tt = t; tt < T[i] - 1; tt++) {
                        // the element that is in R' with the following in null state
                        if (resPaths[seqIndex][tt] == 6 &&
                                                resPaths[seqIndex][tt+1] == 7) {
                            outEnd[k] = tt + 1;
                            break;
                        }
                        // the last element in the sequence is part of the attC site
                        if (tt == (T[i] - 2) && resPaths[seqIndex][T[i]-1] != 7) {
                            outEnd[k] = T[i];
                        }
                    }

                    k++;
                    // continue the search for new attC sites on the remaining
                    // elements in the sequence
                    t = tt + 1;
                } else { // element is not in R''
                    t++;
                }
            }

            seqIndex = 2 * i + 1; // complementary sequence (odd)
            t = 1;

            // check if first element is R'' in an attC site
            if (resPaths[seqIndex][0] == 0) {
                revHits++;
                // increase allocated memory if necessary
                if (k >= size) {
                    size = size * 2;
                    if (!(tmpVec = realloc(outN, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outN = tmpVec;
                    if (!(tmpVec = realloc(outIni, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outIni = tmpVec;
                    if (!(tmpVec = realloc(outEnd, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    outEnd = tmpVec;
                    if (!(tmpVec = realloc(revIndex, size * sizeof(int)))) {
                        printf("Error allocating memory!");
                        return -1;
                    }
                    revIndex = tmpVec;
                }
                // store the sequence number and the first element of attC site
                outN[k] = i + 1;
                outIni[k] = 1;
                revIndex[k] = 1; // reversed is not true

                // find the last element of the attC site and store it
                for (tt = 0; tt < T[i] - 1; tt++) {
                    // the element that is in R' with the following in null state
                    if (resPaths[seqIndex][tt] == 6 &&
                                                resPaths[seqIndex][tt+1] == 7) {
                        outEnd[k] = tt + 1;
                        break;
                    }
                    // the last element in the sequence is part of the attC site
                    if (tt == (T[i] - 2) && resPaths[seqIndex][T[i]-1] != 7) {
                        outEnd[k] = T[i];
                    }
                }
                k++;
                // continue the search for new attC sites on the remaining
                // elements in the sequence
                t = tt + 1;
            }
            // continue searching the sequence for attC sites, by looking for
            // R'' in the remaining elements
            while (t < T[i]) {
                // check if previous element is the null state and current element
                // is R'', if so, new attC site found
                if (resPaths[seqIndex][t-1] == 7 && resPaths[seqIndex][t] == 0) {
                    // increase memory if necessary
                    if (k >= size) {
                        size = size * 2;
                        if (!(tmpVec = realloc(outN, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outN = tmpVec;
                        if (!(tmpVec = realloc(outIni, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outIni = tmpVec;
                        if (!(tmpVec = realloc(outEnd, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        outEnd = tmpVec;
                        if (!(tmpVec = realloc(revIndex, size * sizeof(int)))) {
                            printf("Error allocating memory!");
                            return -1;
                        }
                        revIndex = tmpVec;
                    }
                    revHits++;
                    // store sequence number and first element
                    outN[k] = i + 1;
                    outIni[k] = t + 1;
                    revIndex[k] = 1; // reversed = false
                    // find the last element of the attC site and store it
                    for (tt = t; tt < T[i] - 1; tt++) {
                        // the element that is in R' with the following in null state
                        if (resPaths[seqIndex][tt] == 6 &&
                                                resPaths[seqIndex][tt+1] == 7) {
                            outEnd[k] = tt + 1;
                            break;
                        }
                        // the last element in the sequence is part of the attC site
                        if (tt == (T[i] - 2) && resPaths[seqIndex][T[i]-1] != 7) {
                            outEnd[k] = T[i];
                        }
                    }
                    k++;
                    // continue the search for new attC sites on the remaining
                    // elements in the sequence
                    t = tt + 1;
                } else { // element is not in R''
                    t++;
                }
            }
        }
    } // end search of attC sites



    int n;
    char *currSeq;
    if ((hits + revHits) > 0) { // if at least one attC site was found
        // calculate diffV for each attC site
        double *diffV, *emi, manu_aux, emi_aux, logV8, *tmpDiffV;
        int *seqTransformed;
        diffV = malloc((hits + revHits) * sizeof(double));

        for (int i = 0; i < (hits + revHits); i++) {
            n = outN[i] - 1;

            // find the correct sequence (if original or complementary)
            if (revIndex[i] == 0) {
                currSeq = seqs[n];
            } else {
                currSeq = revSeqs[n];
            }

            // transform the current sequence to 0,1,2,3 instead
            seqTransformed = malloc(T[n] * sizeof(int));
            for (int j = 0; j < T[n]; j++) {
                if (currSeq[j] == 'A' || currSeq[j] == 'a') {
                    seqTransformed[j] = 0;
                }
                if (currSeq[j] == 'C' || currSeq[j] == 'c') {
                    seqTransformed[j] = 1;
                }
                if (currSeq[j] == 'G' || currSeq[j] == 'g') {
                    seqTransformed[j] = 2;
                }
                if (currSeq[j] == 'T' || currSeq[j] == 't') {
                    seqTransformed[j] = 3;
                }
            }

            // calculate the logV score assuming the null state instead of attC site

            // initiate to the logV before the attC site
            if (outIni[i] > 1) {
                logV8 = resLogV[(1 + reverse) * n + revIndex[i]][7 + 8 * (outIni[i] - 2)];
                // add the log probabilities of staying in null state the entire duration
                for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    logV8 += est_emi[7][seqTransformed[t - 1] +
                                            emi_rows[7] * seqTransformed[t]] +
                                            est_trans[7 * 8 + 7];
                }
            } else {
                logV8 = est_ini[7];
                for (t = outIni[i]; t < outEnd[i]; t++) {
                    logV8 += est_emi[7][seqTransformed[t - 1] +
                                            emi_rows[7] * seqTransformed[t]] +
                                            est_trans[7 * 8 + 7];
                }
            }

            // calculate the diffV, the difference between attC site and staying
            // in state 8
            if (outEnd[i] == (T[n])) { //
                diffV[i] =
                resLogV[(1 + reverse) * n + revIndex[i]][resPaths[(1 + reverse)
                                * n + revIndex[i]][(T[n]-1)] + 8 * (T[n]-1)] -
                    logV8 - est_trans[7 * 8 + 7];
            } else {
                diffV[i] =
                resLogV[(1 + reverse) * n + revIndex[i]][6 + 8 * (outEnd[i]-1)] -
                    logV8 - est_trans[7 * 8 + 7];
            }
            free(seqTransformed);
        } // end calc diffV



			// remove overlapping hits on bottom/top 
			// by checking how many hits on the same sequence and comparing their sta with others stop.
			// cases:	1) sta[i] == sto[j]
			//				2) 												sto[i] == sta[j]
			// 				3) sta[i] == sto[j] 	&& 	sto[i] == sta[j]
			// keep the one with higher diffV
			int m = 0;
			//int M = hits + revHits;
			int mm = m, MM, i, j, anyRev = 0;
	
			// for all hits, skipping those on the same seq
			// below "-1" to assure we check the penultimate, which will in turn check the last one in mm.
			while (m < (hits + revHits - 1)) {
				mm = m + 1;
				MM = 0;

				anyRev = 0;
				if (revIndex[m] == 1 ) {
					anyRev = 1;
				}

				// if hits are on the same sequence
				while ( (mm < (hits + revHits)) && (headers[outN[m] - 1] == headers[outN[mm] - 1]) ) {
					if (revIndex[mm] == 1) {
						anyRev = 1;
					}
					MM = MM + 1;
					mm = mm + 1;
				}
				mm = mm -1;

				// if there is more than one hit on the sequence
				if ((MM != 0) && (anyRev == 1) ) {
					// checking the if hits overlap
					//for (i = m	 ; i <  m + MM; i++) {
					//for (j = i +1; j <= m + MM; j++) {
					i = m;
					while (i < m + MM) {
						n = outN[i] - 1;
						j = i + 1;
						while (j <= m + MM) {
							if ((revIndex[i] == 0 && revIndex[j] == 1) || (revIndex[i] == 1 && revIndex[j] == 0)) {
							// if one of the edges match OR if they overlap inexactly 
								if (outIni[i] == (T[n] - outEnd[j] +1) || (T[n] - outIni[i]+1) == outEnd[j] || outEnd[i] == (T[n] - outIni[j] + 1) || (T[n] - outEnd[i] +1) == outIni[j]			|| (outIni[i] > (T[n] - outIni[j] +1) && outIni[i] < (T[n] - outEnd[j] +1)) || (outEnd[i] > (T[n] - outIni[j] +1) && outEnd[i] < (T[n] - outEnd[j] +1))      || ((T[n] - outIni[j] +1) > outIni[i] && (T[n] - outIni[j] +1) < outEnd[i])				|| ((T[n] - outEnd[j] +1) > outIni[i] && (T[n] - outEnd[j] +1) < outEnd[i])				|| (outIni[j] > (T[n] - outEnd[i] +1)  && outIni[j] < (T[n] - outIni[i] +1))				|| (outEnd[j] > (T[n] - outEnd[i] +1)  && outEnd[j] < (T[n] - outIni[i] +1))				||
	((T[n] - outIni[i] +1) > outIni[j] && (T[n] - outIni[i] +1) < outEnd[j])				|| ((T[n] - outEnd[i] +1) > outIni[j] && (T[n] - outEnd[i] +1) < outEnd[j])				) {
									if (diffV[i] >= diffV[j]) {
										//remove j
										if (revIndex[j] == 0) {
											hits = hits -1;								
										} else {
											revHits = revHits -1;
										}
										for (k = j; k < hits + revHits; k++) {
											revIndex[k] = revIndex[k+1];
											outN[k] = outN[k+1];
											outIni[k] = outIni[k+1];
											outEnd[k] = outEnd[k+1];
											diffV[k] = diffV[k+1];
										}
										if (!(tmpVec = realloc(revIndex, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            revIndex = tmpVec;
										if (!(tmpVec = realloc(outN, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            outN = tmpVec;
										if (!(tmpVec = realloc(outIni, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            outIni = tmpVec;
										if (!(tmpVec = realloc(outEnd, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
										if (!(tmpDiffV = realloc(diffV, (hits + revHits) * sizeof(double)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
										diffV = tmpDiffV;
				            outEnd = tmpVec;
										j = j -1;
										MM = MM -1;
									} else {
										// remove i
										if (revIndex[i] == 0) {
											hits = hits -1;								
										} else {
											revHits = revHits -1;
										}
										for (k = i; k < hits + revHits; k++) {
											revIndex[k] = revIndex[k+1];
											outN[k] = outN[k+1];
											outIni[k] = outIni[k+1];
											outEnd[k] = outEnd[k+1];
											diffV[k] = diffV[k+1];
										}
										if (!(tmpVec = realloc(revIndex, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            revIndex = tmpVec;
										if (!(tmpVec = realloc(outN, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            outN = tmpVec;
										if (!(tmpVec = realloc(outIni, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            outIni = tmpVec;
										if (!(tmpVec = realloc(outEnd, (hits + revHits) * sizeof(int)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
				            outEnd = tmpVec;
										if (!(tmpDiffV = realloc(diffV, (hits + revHits) * sizeof(double)))) {
				              printf("Error allocating memory!");
				              return -1;
				            }
										diffV = tmpDiffV;
										i = i -1;
										MM = MM -1;
									}
								}
							}
							j++;
						}
						i++;
					}
				}
				m = m + MM + 1;
		}
		// end remove overlapping hits

        // find end coordinate of each state of each attC site
        int nextState;
        int *stateCoords;
        stateCoords = malloc((S-1) * (hits + revHits) * sizeof(int));

        // for each attC site, loop through the elements and store the current
        // position whenever next state is different from the current state
        for (int i = 0; i < (hits + revHits); i++) {
            n = outN[i] - 1;
            nextState = 1;
            stateCoords[i * (S-1)] = outIni[i];

            if (revIndex[i] == 0) {
                for (int t = outIni[i]; t < outEnd[i]; t++) {
                    if (resPaths[(1 + reverse) * n][t] == nextState) {
                        stateCoords[i * (S - 1) + nextState] = t + 1;
                        nextState++;
                        if (nextState > 6) {
                            break;
                        }
                    }
                }
            } else {
                for (int t = outIni[i]; t < outEnd[i]; t++) {
                    if (resPaths[(1 + reverse) * n + 1][t] == nextState) {
                        stateCoords[i * (S - 1) + nextState] = t + 1;
                        nextState++;
                        if (nextState > 6) {
                            break;
                        }
                    }
                }
            }
        }

        // open files for output
        outfile = fopen(output, "a");
        outTmp = fopen(DIR "/output/out.tmp", "a");
        outTmpFasta = fopen("outHattCI.fasta", "a");

		  char line[150]; 

        // print sequence number, sequence header, first and last element,
        // state coordinates, diffV and indicator of top or bottom strand to file
        for (int i = 0; i < (hits + revHits); i++) {
            n = outN[i] - 1;
            if (revIndex[i] == 0) {
								
                fprintf(outfile,
                "%5i\t%30s\t%5i\t%5i\t%10.8lf\t%5i\t%8i\t%5i\t%5i\t%5i\t%7i\t%5i\t%i\n",
                prevRead + outN[i], headers[outN[i] - 1], outIni[i], outEnd[i],
                diffV[i], stateCoords[i * (S-1)], stateCoords[i * (S-1) + 1],
                stateCoords[i * (S-1) + 2], stateCoords[i * (S-1) + 3],
                stateCoords[i * (S-1) + 4], stateCoords[i * (S-1) + 5],
                stateCoords[i * (S-1) + 6], revIndex[i]);
            } else {
                fprintf(outfile,
                "%5i\t%35s\t%5i\t%5i\t%10.8lf\t%5i\t%8i\t%5i\t%5i\t%5i\t%7i\t%5i\t%i\n",
                prevRead + outN[i], headers[outN[i] - 1], T[n] - outIni[i] + 1, T[n] - outEnd[i] + 1,
                diffV[i], T[n] - stateCoords[i * (S-1)] +1, T[n] - stateCoords[i * (S-1) + 1] +1,
                T[n] - stateCoords[i * (S-1) + 2] +1, T[n] - stateCoords[i * (S-1) + 3] +1,
                T[n] - stateCoords[i * (S-1) + 4] +1, T[n] - stateCoords[i * (S-1) + 5] +1,
                T[n] - stateCoords[i * (S-1) + 6] +1, revIndex[i]);
            }
        }

        // print further information about each attC site to file
        char states[8] = {'R', 's', 'L', '.', 'L', 's', 'R', 'n'};
        for (i = 0; i < (hits + revHits); i++) {
            n = outN[i] - 1;
            if (revIndex[i] == 0) {
                // printing > line in fasta
		    fprintf(outTmpFasta, ">%s_%i_%i Vscore: %lf\n",
                    headers[n], outIni[i], outEnd[i],diffV[i]);
					 // printing attC site
					 for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    fprintf(outTmp, "%c", seqs[n][t]);
						  fprintf(outTmpFasta, "%c", seqs[n][t]);
                }
                fprintf(outTmp, "\n");
                fprintf(outTmpFasta, "\n\n");
					 // printing attC site annotation (path)
                for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    fprintf(outTmp, "%c", states[resPaths[(1 + reverse) * n][t]]);
						  //fprintf(outTmpFasta, "%c", states[resPaths[(1 + reverse) * n][t]]);
                }
                fprintf(outTmp,
                    "\nID: %s\nstart position: %i\nend position: %i\nVscore: %lf\nhit #: %i\n\n",
                    headers[n], outIni[i], outEnd[i],diffV[i], *totHits + *totHitsRev + i+1);
            } else {
                // printing > line in fasta
		    fprintf(outTmpFasta, ">%s_%i_%i Vscore: %lf REVERSED\n",
                    headers[n], T[n]-outIni[i]+1, T[n]-outEnd[i]+1,diffV[i]);
   	        // writing the reversed in the FASTA file, correct so that it appears in top strand form!!
                //for (t = T[n]-outEnd[i]+1-1; t < T[n]-outIni[i]+1; t++) {
		//				  fprintf(outTmpFasta, "%c", seqs[n][t]);
                //}
                for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    fprintf(outTmpFasta, "%c", revSeqs[n][t]);
                }
	        // writing reversed in the outTmp file
                for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    fprintf(outTmp, "%c", seqs[n][T[n] - t - 1]);
                }
                fprintf(outTmp, "\n");
                for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    fprintf(outTmp, "%c", revSeqs[n][t]);
                }
                fprintf(outTmp, "\n");
	        fprintf(outTmpFasta, "\n\n");
                for (t = outIni[i] - 1; t < outEnd[i]; t++) {
                    fprintf(outTmp, "%c", states[resPaths[(1 + reverse) * n + 1][t]]);
                }
                fprintf(outTmp,
                    "\nID: %s\nstart position: %i\n end position: %i\nVscore: %lf\nhit #: %i\nReversed\n\n",
                    headers[n], T[n]-outIni[i]+1, T[n]-outEnd[i]+1, diffV[i], *totHits + *totHitsRev + i+1);
            }
        }
        // close output files
        fclose(outfile);
        fclose(outTmp);
			  fclose(outTmpFasta);

        // add number of hits to counters of total number of hits
        *totHits += hits;
        *totHitsRev += revHits;

        // free variables
        free(diffV);
        free(stateCoords);
    }

    // free variables
    free(revIndex);
    free(outN);
    free(outIni);
    free(outEnd);
    return 0;
}
