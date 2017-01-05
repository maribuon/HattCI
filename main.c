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

Main function for HattCI.

Input arguments:
	- input file: path to fasta file
	- output file: path to output. output file is created if it does not already
		exist. if it exists, the file is overwritten.
	- flags: -b (check both the original sequences and the reversed complementary)
			 -s x (specify chunk size to read and process, x is an int)
			 -t x (specify number of threads to run with)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>

int readFasta(char *, char *, int, double *, double **, double *, double **,
				int *, int *, int, int);
void paraRead(double **, double **, double ***, double ***, int *, int *,
				int *, int);

int main(int argc, char **argv) {
	const int S = 8;
	char **headers, **seqs, *tmp;
	double *est_ini, *est_trans, **f, **est_emi;
	int Dmin[8] = {7, 3, 8, 6, 7, 3, 7, 1},
	 	Dmax[8] = {7, 7, 8, 120, 7, 12, 7, 1},
		Dmax_cum[8], nseq, *T, i, j, skip, status, reverse = 0, maxReads = 1000,
		index, c, base = 0, nthreads = 1;
	char *endptr;

	// check options for running both directions (-b) and number of sequences
	// to read at a time (-s)
	while ((c = getopt (argc, argv, "bs:t:help:version")) != -1) {
    	switch (c) {
      		case 'h' :
        		printf(	"# HattCI :: Hidden markov model for Attc site Identification. \n");
        		printf(	"# HattCI 1.0b (September, 2016) \n");
						printf(	"# Freely distributed under the GNU General Public License (GPLv3). \n\n" );
						printf(	"# - - - - - - - - - - - - - - - - - - - - - - - - - - - \n\n");
						printf(	"Usage: hattci.out <fastaFile> [outFile] [-options] \n\n");
						printf(	"  outFile	: Path to output. Note that if output.txt already exists, it is overwritten. Otherwise, it is created. If no output file is given, the output is written to a standard file \n\n");
						printf(	"  -b		: Processes both the ordinary sequences and the complementary sequence.\n");
						printf(	"  -s [x]	: HattCI reads x sequences at a time and processes them before reading the next x sequences, in order to avoid overextending RAM. Default is 1000 sequences.\n");
						printf(	"  -t [x]	: HattCI may run a large part of the computations in parallel, i.e. let different threads process a set of sequences, which in turn gives a reduced computation time. The parallelization works best when processing larger chunks of sequences at a time. Default number of threads are 1. To specify the number of threads, type for instance. \n\n");
        		return 1;
      		case 'v' :
        		printf(	"# HattCI :: Hidden markov model for Attc site Identification. \n");
        		printf(	"# HattCI 1.0b (September, 2016) \n");
						printf(	"# Freely distributed under the GNU General Public License (GPLv3). \n\n" );
						printf(	"# - - - - - - - - - - - - - - - - - - - - - - - - - - - \n\n");
						printf(	"Cite as: \n\n");
						printf(	"  HattCI: Fast and Accurate attC site Identification Using Hidden Markov Models \n");
						printf(	"  Pereira, Mariana Buongermino; Wallroth, Mikael; Kristiansson, Erik; and Axelson-Fisk, Marina. Journal of Computational Biology. November 2016, 23(11): 891-902. doi:10.1089/cmb.2016.0024.\n");
        		return 1;
      		case 'b' :
        		reverse = 1;
        		break;
      		case 's' :
        		maxReads = strtol(optarg, &endptr, 0);
        		break;
			case 't' :
				nthreads = strtol(optarg, &endptr, 0);
				break;
      		case '?' :
        		if (optopt == 's' || optopt == 't')
          			fprintf(stderr, "Option -%c requires an argument.\n",
							optopt);
        		else if (isprint (optopt))
          			fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        		else
          			fprintf(stderr, "Unknown option character `\\x%x'.\n",
							optopt);
        			return 1;
      		default:
        		abort ();
      	}
	}

	omp_set_num_threads(nthreads);

	// Check argument for input file
	if ((argc - optind) < 1) {
		printf("Input file must be supplied! \n");
		return -1;
	}
	// Check argument for output file
	if ((argc - optind) != 2) {
		printf("No output file was selected, writing results to results.txt\n");
		argv[optind + 1] = "results.txt";
	}

	// initiate pointers to null
	headers = NULL;
	seqs = NULL;
	tmp = NULL;
	est_ini = NULL;
	est_trans = NULL;
	f = NULL;
	est_emi = NULL;

	// load parameters
	Dmax_cum[0] = Dmax[0];
	for (i = 1; i < S; i++) {
		Dmax_cum[i] = Dmax_cum[i-1] + Dmax[i];
	}
	paraRead(&est_ini, &est_trans, &f, &est_emi, Dmin, Dmax, Dmax_cum, S);

 	// call readFasta
	if ((status = readFasta(argv[optind], argv[optind + 1], S, est_ini, est_emi,
		est_trans, f, Dmin, Dmax, reverse, maxReads)) != 0) {
			printf("Error in readFasta!\n");
	}

	// free parameters
	free(est_ini);
	free(est_trans);
	for (i = 0; i < S; i++) {
		free(f[i]);
		free(est_emi[i]);
	}
	free(f);
	free(est_emi);

	return(0); // success
}
