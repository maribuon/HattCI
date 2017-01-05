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

Reads fasta files and runs hmm (wrapper for viterbi and calcResults) on a
specified number of sequences at a time. Returns when all sequences in fasta
file is read and processed.

Parameters:
input 		(input)     char array, fasta file to be processed
output 		(input)     char array, output file to print results to
S			(input)     int, number of states
est_ini		(input)		double array of length S, initial probabilities
est_emi		(input)		double array of arrays, emission probabilities
est_trans	(input)	    double array, transition probabilities
f			(input)		double array of arrays,
Dmin		(input)		double array, min length of each state
Dmax		(input)		double array, max length of each state
reverse     (input)     int, 1 if both directions should be checked, else 0
maxReads    (input)     int, maximum number of sequences to be procssed at a time
*/

#define _POSIX_SOURCE
#define _POSIX_C_SOURCE 200809L
#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // for timing



int hmm(char **, int, double *, double **, double *, double **, int *, int *,
            int *, int, char **, int, int *, int *, int, char *);

int readFasta(char input[], char output[], int S, double *est_ini,
    double **est_emi, double *est_trans, double **f, int *Dmin, int *Dmax,
    int reverse, int maxReads) {

		// random number generator seed
		srand48(time(NULL));

    FILE *infile;
    ssize_t bytes_read;
    size_t nbytes = 0;
    char *lineptr = NULL, **seqs, **headers, *tmpVec;
    int rows = 0, skip = 0, i, j, *T, seqrows, prevRead = 0, totHits = 0,
        totHitsRev = 0, notProcessed = 0, seqsProcessed = 0, status;
    FILE *outfile;
    FILE *outTmp;
    long int mbases = 0;
		float pn;




    infile = fopen(input,"r");
    outfile = fopen(output, "w");

    // write header of output file
		// diffV now called Vscore
    fprintf(outfile,
        "%5s\t%35s\t%5s\t%5s\t%10s\t%5s\t%8s\t%5s\t%5s\t%5s\t%7s\t%5s\n", "hit", "ID", "start",
        "end", "Vscore", "R\'\'", "spacer\'\'", "L\'\'", "loop", "L\'", "spacer\'", "R\'");
    fclose(outfile);

    // allocate memory to store headers, seqs and sequence lengths
    headers = malloc(maxReads * sizeof(char *));
    seqs = malloc(maxReads * sizeof(char *));
    T = malloc(maxReads * sizeof(int));

    // ensure allocations were successful
    if (!headers || !seqs || !T) {
        printf("Error allocating memory!\n");
        return -1;
    }

    int nseq = -1; // counter for sequences to be processed in each chunk
    while ((bytes_read = getline(&lineptr, &nbytes, infile)) != -1) {
        // if read line is a header, new sequence found to be read and processed
        if (lineptr[0] == '>') {
            nseq++;
            seqsProcessed++;

            // if we already have filled the current chunk to process, then
            // process the sequences
            if (nseq >= maxReads) {
                for (i = 0; i < nseq ; i++) {
                    skip = 0;

                    // replace all 'n' (error in sequence) by random {A,C,G,T} 
										// and remove all new line characters from the sequences
                    for (j = 0; j < T[i]; j++) {
                        if (seqs[i][j] == 'A' || seqs[i][j] == 'a' ||
                            seqs[i][j] == 'C' || seqs[i][j] == 'c' ||
                            seqs[i][j] == 'G' || seqs[i][j] == 'g' ||
                            seqs[i][j] == 'T' || seqs[i][j] == 't') {
                            seqs[i][j-skip] = seqs[i][j];
                        } else if (seqs[i][j] == 'N') {
														pn = drand48();
														if (pn <= 0.25) {
                            	seqs[i][j-skip] = 'A';
														}
														else if (pn <= .5) {
                            	seqs[i][j-skip] = 'C';
														}
														else if (pn <= .75) {
                            	seqs[i][j-skip] = 'G';
														} 
														else {
                            	seqs[i][j-skip] = 'T';
														}
                        } else if (seqs[i][j] == 'n') {
													pn = drand48();
														if (pn <= 0.25) {
                            	seqs[i][j-skip] = 'a';
														}
														else if (pn <= .5) {
                            	seqs[i][j-skip] = 'c';
														}
														else if (pn <= .75) {
                            	seqs[i][j-skip] = 'g';
														} 
														else {
                            	seqs[i][j-skip] = 't';
														}
                        } else {
                            skip++;
                        }
                    }
                    T[i] = T[i] - skip; // the "true" length of the sequence


                    // reallocate to the "true" size, check if successful
                    if (T[i] > 0) {
                        if (!(tmpVec = realloc(seqs[i], T[i]))) {
                            printf("Error reallocating memory for sequence (reducing), %i!\n", T[i]);
                            return -1;
                        }
                    } else {
                        T[i] = 1;
                        if (!(tmpVec = realloc(seqs[i], T[i]))) {
                            printf("Error reallocating memory for sequence (reducing), %i!\n", T[i]);
                            return -1;
                        }
                    }
                    seqs[i] = tmpVec;


                    j = 0;
                    // store only the info in the header up to the first
                    // white space
                    while(headers[i][j+1] != ' ' && headers[i][j+1] != '\n') {
                        j++;
                        headers[i][j-1] = headers[i][j];
                    }
						  
                    // reallocate to the new size of the header and including
                    // a null character at the end
                    if (!(tmpVec = realloc(headers[i],
                                                (j + 1) * sizeof(char)))) {
                        printf("Error reallocating memory for header (reducing)!\n");
                        return -1;
                    }
                    headers[i] = tmpVec;
                    headers[i][j] = '\0';
                    mbases += T[i];
                }

                // call on hmm to process the sequences
                if ((status = hmm(seqs, S, est_ini, est_emi, est_trans, f, Dmin,
                    Dmax, T, nseq, headers, prevRead, &totHits, &totHitsRev,
                    reverse, output)) != 0) {
                    printf("Error in hmm\n");
                    return -1;
                }

                // reset the list of sequences and headers (free memory)
                for (i = 0; i < nseq; i++) {
                    free(headers[i]);
                    free(seqs[i]);
                }
                prevRead += nseq; // store the total number of sequences read this far
                nseq = 0; // number of sequences read in the new chunk
            }
            // store the read header
            headers[nseq] = malloc(bytes_read * sizeof(char));
            memcpy(headers[nseq], lineptr, bytes_read);
            seqrows = 0;
        } else if (lineptr[0] != '\n') { // read line is (part of) a sequence
            if (seqrows == 0) { // first line of a sequence
                if (!(seqs[nseq] = malloc(bytes_read * sizeof(char)))) {
                    printf("Error allocating memory for sequence!");
                    return -1;
                }
                memcpy(seqs[nseq], lineptr, bytes_read);
                T[nseq] = bytes_read;
                seqrows++;
            } else { // append to current sequence
                if (!(tmpVec = realloc(seqs[nseq],
                                (bytes_read + T[nseq] + 1) * sizeof(char)))) {
                    printf("Error reallocating memory for sequence (increasing)!");
                    return -1;
                }
                seqs[nseq] = tmpVec;
                seqs[nseq][T[nseq]] = '\0';
                strcat(seqs[nseq], lineptr);
                T[nseq] += bytes_read;
                seqrows++;
            }
        } // otherwise we have an empty line, do nothing

        // prepare for reading next line
        free(lineptr);
        lineptr = NULL;
        nbytes = 0;
    } // end while (end of file)

    free(lineptr);
    lineptr = NULL;
    nbytes = 0;
    fclose(infile);

    nseq++; // set nseq to number of seqs instead of index (0, 1, ..., n-1)
    for (i = 0; i < nseq ; i++) {
        skip = 0;

        // remove all 'n' (error in sequence) and new line characters from the sequences
        for (j = 0; j < T[i]; j++) {
            if (seqs[i][j] == 'A' || seqs[i][j] == 'a' ||
                seqs[i][j] == 'C' || seqs[i][j] == 'c' ||
                seqs[i][j] == 'G' || seqs[i][j] == 'g' ||
                seqs[i][j] == 'T' || seqs[i][j] == 't') {
                seqs[i][j-skip] = seqs[i][j];
            } else if (seqs[i][j] == 'N') {
								pn = drand48();
								if (pn <= 0.25) {
		            	seqs[i][j-skip] = 'A';
								}
								else if (pn <= .5) {
		            	seqs[i][j-skip] = 'C';
								}
								else if (pn <= .75) {
		            	seqs[i][j-skip] = 'G';
								} 
								else {
		            	seqs[i][j-skip] = 'T';
								}
		        } else if (seqs[i][j] == 'n') {
							pn = drand48();
								if (pn <= 0.25) {
		            	seqs[i][j-skip] = 'a';
								}
								else if (pn <= .5) {
		            	seqs[i][j-skip] = 'c';
								}
								else if (pn <= .75) {
		            	seqs[i][j-skip] = 'g';
								} 
								else {
		            	seqs[i][j-skip] = 't';
								}
            } else {
                skip++;
         		}
        }
        T[i] = T[i] - skip; // the "true" length of the sequence


        if (T[i] > 0) {
            if (!(tmpVec = realloc(seqs[i], T[i]))) {
                printf("Error reallocating memory for sequence (reducing), %i!\n", T[i]);
                return -1;
            }
        } else {
            T[i] = 1;
            if (!(tmpVec = realloc(seqs[i], T[i]))) {
                printf("Error reallocating memory for sequence (reducing), %i!\n", T[i]);
                return -1;
            }
        }
        seqs[i] = tmpVec;

        j = 0;
        // store only the info in the header up to the first white space
        while(headers[i][j+1] != ' ' && headers[i][j+1] != '\n') {
					j++;
            	headers[i][j-1] = headers[i][j];
        }
		  

        // reallocate to the new size of the header and including a null
        // character at the end
        if (!(tmpVec = realloc(headers[i], (j + 1) * sizeof(char)))) {
            printf("Error allocating memory!\n");
            return -1;
        }

        headers[i] = tmpVec;
        headers[i][j] = '\0';
        mbases += T[i];

    }

    // call on hmm to process the sequences
    if ((status = hmm(seqs, S, est_ini, est_emi, est_trans, f, Dmin, Dmax, T,
        nseq, headers, prevRead, &totHits, &totHitsRev, reverse, output)) != 0) {
            printf("Error in hmm\n");
            return -1;
        }

    // free memory
    for (i = 0; i < nseq; i++) {
		free(headers[i]);
		free(seqs[i]);
	}
    free(headers);
    free(seqs);
    free(T);

    // write final output (append outTmp to outfile)
    if (totHits + totHitsRev > 0) {
        outfile = fopen(output, "a");
        outTmp = fopen(DIR "/output/out.tmp", "r");

        fprintf(outfile, "--------------------------------------------\n");
        while ((bytes_read = getline(&lineptr, &nbytes, outTmp)) != -1) {
            fprintf(outfile, "%s", lineptr);
            free(lineptr);
            lineptr = NULL;
            nbytes = 0;
        }
        free(lineptr);
        lineptr = NULL;
        nbytes = 0;
        fclose(outTmp);
        fclose(outfile);

        if (remove(DIR "/output/out.tmp") != 0) {
            printf("Could not delete out.tmp\n");
        }


    }
    // print summary of processed sequences to standard output
    printf("Mbases processed: %lf\nSequences read: %i\nSequences processed: %i\n",
            mbases * 1.0 / 1000000, seqsProcessed, prevRead + nseq);
    if (reverse == 0) {
        printf("Hits: %i\n", totHits);
    } else {
        printf("Hits (top strand): %i\nHits (bottom strand): %i\n", totHits,
                totHitsRev);
    }

    return 0; // readfasta successful
}
