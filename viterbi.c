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

Runs the Viterbi algorithm on a sequence.
The states of the algorithm is
- initiation (initiates the starting conditions)
- recursion (calculate the likelihood of each possible path)
- termination (terminates the path by a transition to the null state from the last element in the sequence)
- find path (moves backwards to find the most probable path through the sequence)

Parameters:
path 		(input/output)	int array of length (seqLen + 1)
logV 		(input/output)	double matrix, size S * (seqLen + 1)
seqs 		(input)			char array of length seqLen
S			(input)			int, number of states
est_ini		(input)			double array of length S, initial probabilities
est_emi		(input)			double array of arrays, emission probabilities
est_trans	(input)			double array, transition probabilities
f			(input)			double array of arrays,
Dmin		(input)			double array, min length of each state
Dmax		(input)			double array, max length of each state
seqLen		(input)			int, length of the sequence

STATES:
0 - R''
1 - spacer''
2 - L''
3 - loop (gHMM)
4 - L'
5 - spacer'
6 - R'
7 - silent (null state)
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void viterbi(int *path, double *logV, char *seqs, int S, double *est_ini,
	double **est_emi, double *est_trans, double **f, int *Dmin, int *Dmax,
	int seqLen) {

	int **psi, **phi, imax, emi_rows[8] = {1, 4, 4, 16, 4, 4, 1, 4};
	double emi_cum, tomax[2] = {0, 0}, *tomax2, tomax_d;
	int curr_state, prev_state, posmax, dmax, length, kk, Dmax_local, d_aux, d,
		kstop, *seq, t;

	// allocate memory for psi (prev state), phi (state duration) and
	// int representation of sequence
	psi = malloc(S * sizeof(int *));
	phi = malloc(S * sizeof(int *));
	seq = malloc(seqLen * sizeof(int));
	for (int i = 0; i < S; i++) {
		psi[i] = malloc((seqLen + 1) * sizeof(int));
		phi[i] = malloc((seqLen + 1) * sizeof(int));
	}

	// initiate logV, phi and psi
	for (int i = 0; i < S - 1; i++) {
		logV[i] = -INFINITY;
		phi[i][0] = -1;
		psi[i][0] = -1;
	}

	// assume the sequence may start in null state
	logV[S-1] = (est_ini[S-1]);
	curr_state = 7;

	// transform the first nucleotides to int repr.
	for (int i = 0; i < 7; i++) {
		if (seqs[i] == 'A' || seqs[i] == 'a') {
			seq[i] = 1;
		} else if (seqs[i] == 'C' || seqs[i] == 'c') {
			seq[i] = 2;
		} else if (seqs[i] == 'G' || seqs[i] == 'g') {
			seq[i] = 3;
		} else {
			seq[i] = 4;
		}
	}

	// assume null state in 7 first elements
	for (t = 1; t < 7; t++) {
		tomax[0] = logV[6 + 8 * (t-1)] + (est_trans[6 * 8 + curr_state]);
		tomax[1] = logV[7 + 8 * (t-1)] + (est_trans[7 * 8 + curr_state]);
		if (tomax[0] >= tomax[1]) {
			imax = 6;
		} else {
			imax = 7;
		}

		logV[curr_state + 8 * t] = tomax[imax-6] +
								(est_emi[curr_state][(seq[t-1]-1) +
								emi_rows[curr_state] * (seq[t]-1)]);
		psi[curr_state][t] = imax;
		phi[curr_state][t] = 1;
		for (int i = 0; i < S - 1; i++) {
			logV[i + 8 * t] = -INFINITY;
		}
	}

	// assume sequence may begin in R''
	t = 6;
	curr_state = 0;
	double logemi;

	// check necessary condition xxxxAAC for R'' state
	if (seq[t-2] != 1 || seq[t-1] != 1 || seq[t] != 2) {
		logemi = -INFINITY;
	} else {
		logemi = (est_emi[curr_state][64*(seq[t-6]-1) + 16*(seq[t-5]-1) +
							4*(seq[t-4]-1) + seq[t-3] - 1]);
	}

	// calculate logV and assign psi and phi, assuming sequence begins in R''
	// assign initiating probabilities from state 0-6 to state 0 (note log prob)
	double sum = 0;
	for (int i = 0; i < S - 1; i++) {
		sum += exp(est_ini[i]);
	}
	logV[curr_state + 8 * t] = log(sum) + logemi;
	psi[curr_state][t] = -1;
	phi[curr_state][t] = 7;

	/************ BEGIN RECURSION *************/

	// calculate minimum first element for each state
	int Dmin_cum[8];
	Dmin_cum[0] = Dmin[0] - 1;
	int Dmin_sum = Dmin[0] - 1;
	for (int i = 1; i < 7; i++) {
		Dmin_sum += Dmin[i];
		Dmin_cum[i] = Dmin_sum;
	}
	Dmin_cum[7] = Dmin_cum[6] + Dmin[7];

	/*
	iterate over each element up to the minimum length of an attC site
	for each element, calculate logV assuming the current element is the last
	element of each state, and assign previous state and duration (psi, phi)
	*/
	for (t = 7; t < Dmin_sum; t++) {
		// transform current element to integer representation
		if (seqs[t] == 'A' || seqs[t] == 'a') {
			seq[t] = 1;
		} else if (seqs[t] == 'C' || seqs[t] == 'c') {
			seq[t] = 2;
		} else if (seqs[t] == 'G' || seqs[t] == 'g') {
			seq[t] = 3;
		} else {
			seq[t] = 4;
		}

		// state 0 (R'') : d = 7, prev state = 7
		curr_state = 0;
		prev_state = 7;
		d = 7;
		logemi = -INFINITY; // base case

		// check necessary condition xxxxAAC and calculate emission probability
		if (seq[t-2] == 1 && seq[t-1] == 1 && seq[t] == 2) {
			logemi = est_emi[curr_state][64*(seq[t-6]-1) + 16*(seq[t-5]-1) +
							4*(seq[t-4]-1) + seq[t-3] - 1];
		}

		// calculate and store logV, psi and phi
		logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-7)] +
							  (est_trans[prev_state * 8 + curr_state]) + logemi;
		psi[curr_state][t] = prev_state;
		phi[curr_state][t] = d;

		// state 1 (SPACER''): d = 4-6, prev state = 0
		curr_state = 1;

		// check if an attC site could fit
		if (t < Dmin_cum[curr_state]) {
			// attC site could not fit
			logV[curr_state + 8 * t] = -INFINITY;
			psi[curr_state][1] = -1;
			phi[curr_state][1] = 0;
		} else {
			// attC site could fit, assign maximum possible state duration
			if (t < Dmax[curr_state]) {
				Dmax_local = t;
			} else {
				Dmax_local = Dmax[curr_state];
			}

			// initiate logV to base case for each possible duration
			prev_state = 0;
			length = Dmax_local - Dmin[curr_state] + 1;
			tomax2 = malloc(length * sizeof(double));
			for (int i = 0; i < length; i++) {
				tomax2[i] = -INFINITY; // base case
			}

			// calculate logV for each possible duration

			// calculate logV for minimum duration
			d_aux = 0;
			emi_cum = 0;
			d = Dmin[curr_state] - 1;
			for (int k = t; k > t - Dmin[curr_state]; k--) {
				emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
							emi_rows[curr_state] * (seq[k]-1)]);
			}
			tomax_d = emi_cum + (f[curr_state][d]);
			tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
						   logV[prev_state + 8 * (t-d-1)] + tomax_d;
			d_aux++;
			posmax = 0;

			if (t < Dmax[curr_state] + 1) {
				kstop = 0;
			} else {
				kstop = t - Dmax[curr_state];
			}

			// calculate logV for the remaining possible durations
			for (int k = t - Dmin[curr_state]; k > kstop; k--) {
				d++;
				emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
						    emi_rows[curr_state] * (seq[k]-1)]);
				tomax_d = emi_cum + (f[curr_state][d]);
				tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
							   logV[prev_state + 8 * (t-d-1)] + tomax_d;
				if (tomax2[d_aux] > tomax2[posmax]) {
					posmax = d_aux; // store the currently most probable case
				}
				d_aux++;
			}

			// store values of the most probable duration of state
			logV[curr_state + 8 * t] = tomax2[posmax];
			psi[curr_state][t] = prev_state;
			phi[curr_state][t] = posmax + Dmin[curr_state];
			free(tomax2);
		}

		// state 2 (L''): d = 8, prev state = 1
		curr_state = 2;

		// check if attC site could fit
		if (t < Dmin_cum[curr_state]) {
			logV[curr_state + 8 * t] = -INFINITY;
			psi[curr_state][1] = -1;
			phi[curr_state][1] = 0;
		} else {
			// attC site could fit, calculate the logV assuming the state is L''
			prev_state = 1;
			d = 8;
			logemi = 0;
			kk = 7;
			// calculate emission probabilities
			for (int k = 0; k < 8; k++) {
				logemi += (est_emi[curr_state][(seq[t-k]-1) +
				 				emi_rows[curr_state] * kk]);
				kk--;
			}

			// calculate logV and store logV, psi and phi
			logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-8)] +
							(est_trans[prev_state * 8 + curr_state]) + logemi;
			psi[curr_state][t] = prev_state;
			phi[curr_state][t] = d;
		}

		// state 3 (loop, gHMM): d = 20-130, prev state = 2
		curr_state = 3;

		// check if attC site could fit
		if (t < Dmin_cum[curr_state]) {
			logV[curr_state + 8 * t] = -INFINITY;
			psi[curr_state][1] = -1;
			phi[curr_state][1] = 0;
		} else {
			// attC site could fit, assign maximum possible duration of state
			prev_state = 2;
			length = Dmax[curr_state] - Dmin[curr_state] + 1;

			// initiate to base case for each possibility
			tomax2 = malloc(length * sizeof(double));
			for (int i = 0; i < length; i++) {
				tomax2[i] = -INFINITY;
			}

			// calculate the logV for each possible duration

			// calculate the emission for the minimum duration of state
			d_aux = 0;
			emi_cum = 0;
			d = Dmin[curr_state] - 1;
			for (int k = t; k > t - Dmin[curr_state]; k--) {
				 emi_cum += est_emi[curr_state][(4*(seq[k-2]-1) +
				 			seq[k-1] - 1) + emi_rows[curr_state] * (seq[k]-1)];
			}

			// calculate logV for the minimum duration of state
			tomax_d = emi_cum + (f[curr_state][d]);
			tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
							logV[prev_state + 8 * (t-d-1)] + tomax_d;
			d_aux++;
			posmax = 0;

			if (t < Dmax[curr_state] + 2) {
				kstop = 1;
			} else {
				kstop = t - Dmax[curr_state];
			}


			// calculate logV for the remaining possible durations
			for (int k = t - Dmin[curr_state]; k > kstop; k--) {
				d++;
				emi_cum += (est_emi[curr_state][(4*(seq[k-2]-1) +
							seq[k-1] - 1) + emi_rows[curr_state] * (seq[k]-1)]);
				tomax_d = emi_cum + (f[curr_state][d]);
				tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
							   logV[prev_state + 8 * (t-d-1)] + tomax_d;

				if (tomax2[d_aux] > tomax2[posmax]) {
					posmax = d_aux; // store the most probable case
				}
				d_aux++;
			}

			// store values of the most probable duration of state
			logV[curr_state + 8 * t] = tomax2[posmax];
			psi[curr_state][t] = prev_state;
			phi[curr_state][t] = posmax + Dmin[curr_state];
			free(tomax2);
		}

		// state 4 (L'): d = 7, prev state = 3
		curr_state = 4;

		// check that attC site could fit
		if (t < Dmin_cum[curr_state]) {
			logV[curr_state + 8 * t] = -INFINITY;
			psi[curr_state][1] = -1;
			phi[curr_state][1] = 0;
		} else {
			// calculate logV for state and store the values
			prev_state = 3;
			d = 7;
			logemi = 0;
			int kk = 6;
			for (int k = 0; k < 7; k++) {
				logemi += (est_emi[curr_state][(seq[t-k]-1) +
							emi_rows[curr_state] * kk]);
				kk--;
			}

			logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-7)] +
							(est_trans[prev_state * 8 + curr_state]) + logemi;
			psi[curr_state][t] = prev_state;
			phi[curr_state][t] = d;
		}

		// state 5 (SPACER'', gHMM): d = 3-12, prev state = 4
		curr_state = 5;

		// check that attC site could fitt
		if (t < Dmin_cum[curr_state]) {
			logV[curr_state + 8 * t] = -INFINITY;
			psi[curr_state][1] = -1;
			phi[curr_state][1] = 0;
		} else {
			// attC site could fit, find the maximum possible duration of the state
			if (t < Dmax[curr_state]) {
				Dmax_local = t;
			} else {
				Dmax_local = Dmax[curr_state];
			}

			// initiate logV to base case for each possible duration
			prev_state = 4;
			length = Dmax_local - Dmin[curr_state] + 1;
			tomax2 = malloc(length * sizeof(double));
			for (int i = 0; i < length; i++) {
				tomax2[i] = -INFINITY;
			}

			// calculate logV for each possible duration

			// calculate logV for minimum duration
			d_aux = 0;
			emi_cum = 0;
			d = Dmin[curr_state] - 1;

			for (int k = t; k > t - Dmin[curr_state]; k--) {
				emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
										emi_rows[curr_state] * (seq[k]-1)]);
			}
			tomax_d = emi_cum + (f[curr_state][d]);
			tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
						   logV[prev_state + 8 * (t-d-1)] + tomax_d;
			d_aux++;
			posmax = 0;

			if (t < Dmax[curr_state] + 1) {
				kstop = 0;
			} else {
				kstop = t - Dmax[curr_state];
			}

			// calculate logV for the remaining possible durations
			for (int k = t - Dmin[curr_state]; k > kstop; k--) {
				d++;
				emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
										 emi_rows[curr_state] * (seq[k]-1)]);
				tomax_d = emi_cum + (f[curr_state][d]);
				tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
							   logV[prev_state + 8 * (t-d-1)] + tomax_d;
				if (tomax2[d_aux] > tomax2[posmax]) {
					posmax = d_aux; // store most probable case
				}
				d_aux++;
			}

			// store values for most probable duration
			logV[curr_state + 8 * t] = tomax2[posmax];
			psi[curr_state][t] = prev_state;
			phi[curr_state][t] = posmax + Dmin[curr_state];
			free(tomax2);
		}

		// state 6 (R') : d = 7, prev state = 5
		curr_state = 6;
		logemi = -INFINITY;

		// check if attC site could fit
		if (t < Dmin_cum[curr_state]) {
			logV[curr_state + 8 * t] = logemi;
			psi[curr_state][1] = -1;
			phi[curr_state][1] = 0;
		} else {
			// attC site could fit, check if necessary condition GTTxxxx holds
			prev_state = 5;
			d = 7;
			if (seq[t-6] == 3 && seq[t-5] == 4 && seq[t-4] == 4) {
				logemi = est_emi[curr_state][64*(seq[t-3]-1) + 16*(seq[t-2]-1) +
							4*(seq[t-1]-1) + seq[t] - 1];
			}

			// calculate and store values for state
			logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-7)] +
							(est_trans[prev_state * 8 + curr_state]) + logemi;
			psi[curr_state][t] = prev_state;
			phi[curr_state][t] = d;

		}

		// state 7 (NULL): d = 1, prev state = {6,7}
		// calculate logV for state assuming that the previous state may be
		// either null state again or R'
		curr_state = 7;
		tomax[0] = logV[6 + 8 * (t-1)] + (est_trans[6 * 8 + curr_state]);
		tomax[1] = logV[7 + 8 * (t-1)] + (est_trans[7 * 8 + curr_state]);

		if (tomax[0] >= tomax[1]) {
			imax = 6;
		} else {
			imax = 7;
		}
		// calculate and store most probable case
		logV[curr_state + 8 * t] = tomax[imax-6] +
								(est_emi[curr_state][(seq[t-1]-1) +
								emi_rows[curr_state] * (seq[t]-1)]);
		psi[curr_state][t] = imax;
		phi[curr_state][t] = 1;
	}

	/*
	continued recursion for all elements larger than the minimum length of attC site
	this implies that checking if attC site could fit is no longer necessary
	*/

    for (t = Dmin_sum; t < seqLen; t++) {
		// transform current element to integer representation
		if (seqs[t] == 'A' || seqs[t] == 'a') {
			seq[t] = 1;
		} else if (seqs[t] == 'C' || seqs[t] == 'c') {
			seq[t] = 2;
		} else if (seqs[t] == 'G' || seqs[t] == 'g') {
			seq[t] = 3;
		} else {
			seq[t] = 4;
		}

        // state 0 (R'') : d = 7, prev state = 7
        curr_state = 0;
        prev_state = 7;
		d = 7;
        logemi = -INFINITY;

		// check necessary condition xxxxAAC and calculate emission probability
		if (seq[t-2] == 1 && seq[t-1] == 1 && seq[t] == 2) {
			logemi = est_emi[curr_state][64*(seq[t-6]-1) + 16*(seq[t-5]-1) +
							4*(seq[t-4]-1) + seq[t-3] - 1];
		}

		// calculate and store logV, psi and phi
        logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-7)] +
        					  (est_trans[prev_state * 8 + curr_state]) + logemi;
        psi[curr_state][t] = prev_state;
        phi[curr_state][t] = d;

		// state 1 (SPACER''): d = 3-7, prev state = 0
		curr_state = 1;
		// assign maximum possible state duration
		if (t < Dmax[curr_state]) {
			Dmax_local = t;
		} else {
			Dmax_local = Dmax[curr_state];
		}

		// initiate logV to base case for each possible duration
		prev_state = 0;
		length = Dmax_local - Dmin[curr_state] + 1;
		tomax2 = malloc(length * sizeof(double));
		for (int i = 0; i < length; i++) {
			tomax2[i] = -INFINITY;
		}

		// calculate logV for each possible duration

		// calculate logV for minimum duration
		d_aux = 0;
		emi_cum = 0;
		d = Dmin[curr_state] - 1;

		for (int k = t; k > t - Dmin[curr_state]; k--) {
			emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
			 			emi_rows[curr_state] * (seq[k]-1)]);
		}

		tomax_d = emi_cum + (f[curr_state][d]);
		tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
					   logV[prev_state + 8 * (t-d-1)] + tomax_d;

		d_aux++;
		posmax = 0;

		if (t < Dmax[curr_state] + 1) {
			kstop = 0;
		} else {
			kstop = t - Dmax[curr_state];
		}

		// calculate logV for the remaining possible durations
		for (int k = t - Dmin[curr_state]; k > kstop; k--) {
			d++;
			emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
						emi_rows[curr_state] * (seq[k]-1)]);
			tomax_d = emi_cum + (f[curr_state][d]);
			tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
						   logV[prev_state + 8 * (t-d-1)] + tomax_d;
			if (tomax2[d_aux] > tomax2[posmax]) {
				posmax = d_aux; // store the currently most probable case
			}

			d_aux++;
		}

		// store values of the most probable duration of state
		logV[curr_state + 8 * t] = tomax2[posmax];
		psi[curr_state][t] = prev_state;
		phi[curr_state][t] = posmax + Dmin[curr_state];
		free(tomax2);

		// state 2 (L''): d = 8, prev state = 1
		curr_state = 2;
		prev_state = 1;
		d = 8;
		logemi = 0;
		kk = 7;

		// calculate the emission probabilities
		for (int k = 0; k < 8; k++) {
			logemi += (est_emi[curr_state][(seq[t-k]-1) +
							emi_rows[curr_state] * kk]);
			kk--;
		}

		// calculate the logV assuming the state is L''
		logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-8)] +
							  (est_trans[prev_state * 8 + curr_state]) + logemi;
		psi[curr_state][t] = prev_state;
		phi[curr_state][t] = d;

		// state 3 (gHMM): d = 20-130, prev state = 2
		curr_state = 3;
		prev_state = 2;

		// assign maximum possible duration of state
		length = Dmax[curr_state] - Dmin[curr_state] + 1;

		// initiate to base case for each possibility
		tomax2 = malloc(length * sizeof(double));
		for (int i = 0; i < length; i++) {
			tomax2[i] = -INFINITY;
		}

		// calculate the logV for each possible duration

		// calculate the emission for the minimum duration of state
		d_aux = 0;
		emi_cum = 0;
		d = Dmin[curr_state] - 1;

		for (int k = t; k > t - Dmin[curr_state]; k--) {
			 emi_cum += (est_emi[curr_state][(4*(seq[k-2]-1) + seq[k-1] - 1) +
									emi_rows[curr_state] * (seq[k]-1)]);
		}

		// calculate logV for the minimum duration of state
		tomax_d = emi_cum + (f[curr_state][d]);
		tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
						logV[prev_state + 8 * (t-d-1)] + tomax_d;
		d_aux++;
		posmax = 0;


		if (t < Dmax[curr_state] + 2) {
			kstop = 1;
		} else {
			kstop = t - Dmax[curr_state];
		}

		// calculate logV for the remaining possible durations
		for (int k = t - Dmin[curr_state]; k > kstop; k--) {
			d++;
			emi_cum += (est_emi[curr_state][(4*(seq[k-2]-1) + seq[k-1] - 1) +
											emi_rows[curr_state] * (seq[k]-1)]);

			tomax_d = emi_cum + (f[curr_state][d]);
			tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
						   logV[prev_state + 8 * (t-d-1)] + tomax_d;

			if (tomax2[d_aux] > tomax2[posmax]) {
				posmax = d_aux; // store the most probable case
			}
			d_aux++;
		}

		// store values of the most probable duration of state
		logV[curr_state + 8 * t] = tomax2[posmax];
		psi[curr_state][t] = prev_state;
		phi[curr_state][t] = posmax + Dmin[curr_state];

		free(tomax2);

		// state 4 (L'): d = 7, prev state = 3
		curr_state = 4;
		prev_state = 3;
		d = 7;

		// calculate logV for state and store the values
		logemi = 0;
		int kk = 6;
		for (int k = 0; k < 7; k++) {
			logemi += (est_emi[curr_state][(seq[t-k]-1) +
								emi_rows[curr_state] * kk]);
			kk--;
		}

		logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-7)] +
							  (est_trans[prev_state * 8 + curr_state]) + logemi;
		psi[curr_state][t] = prev_state;
		phi[curr_state][t] = d;

		// state 5 (SPACER', gHMM): d = 4-8, prev state = 4
		curr_state = 5;
		prev_state = 4;

		//find the maximum possible duration of the state
		if (t < Dmax[curr_state]) {
			Dmax_local = t;
		} else {
			Dmax_local = Dmax[curr_state];
		}

		// initiate logV to base case for each possible duration
		length = Dmax_local - Dmin[curr_state] + 1;
		tomax2 = malloc(length * sizeof(double));
		for (int i = 0; i < length; i++) {
			tomax2[i] = -INFINITY;
		}

		// calculate logV for each possible duration

		// calculate logV for minimum duration
		d_aux = 0;
		emi_cum = 0;
		d = Dmin[curr_state] - 1;

		for (int k = t; k > t - Dmin[curr_state]; k--) {
			emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
						emi_rows[curr_state] * (seq[k]-1)]);
		}
		tomax_d = emi_cum + (f[curr_state][d]);
		tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
					   logV[prev_state + 8 * (t-d-1)] + tomax_d;
		d_aux++;
		posmax = 0;

		if (t < Dmax[curr_state] + 1) {
			kstop = 0;
		} else {
			kstop = t - Dmax[curr_state];
		}

		// calculate logV for the remaining possible durations
		for (int k = t - Dmin[curr_state]; k > kstop; k--) {
			d++;
			emi_cum += (est_emi[curr_state][(seq[k-1]-1) +
							emi_rows[curr_state] * (seq[k]-1)]);
			tomax_d = emi_cum + (f[curr_state][d]);
			tomax2[d_aux] = (est_trans[prev_state * 8 + curr_state]) +
						   logV[prev_state + 8 * (t-d-1)] + tomax_d;
			if (tomax2[d_aux] > tomax2[posmax]) {
				posmax = d_aux; // store most probable case
			}
			d_aux++;
		}

		// store values for most probable duration
		logV[curr_state + 8 * t] = tomax2[posmax];
		psi[curr_state][t] = prev_state;
		phi[curr_state][t] = posmax + Dmin[curr_state];

		free(tomax2);

        // state 6 (R') : d = 7, prev state = 5
        curr_state = 6;
        prev_state = 5;
        logemi = -INFINITY;
		d = 7;

		// check if necessary condition GTTxxxx holds, if so calculate for emission
        if (seq[t-6] == 3 && seq[t-5] == 4 && seq[t-4] == 4) {
			logemi = est_emi[curr_state][64*(seq[t-3]-1) + 16*(seq[t-2]-1) +
							4*(seq[t-1]-1) + seq[t] - 1];
		}

		// calculate and store values for state
        logV[curr_state + 8 * t] = logV[prev_state + 8 * (t-7)] +
                            (est_trans[prev_state * 8 + curr_state]) + logemi;
        psi[curr_state][t] = prev_state;
        phi[curr_state][t] = d;

		// state 7 (NULL): d = 1, prev state = {6,7}
		// calculate logV for state assuming that the previous state may be
		// either null state again or R'
		curr_state = 7;
		tomax[0] = logV[6 + 8 * (t-1)] + (est_trans[6 * 8 + curr_state]);
		tomax[1] = logV[7 + 8 * (t-1)] + (est_trans[7 * 8 + curr_state]);

		if (tomax[0] >= tomax[1]) {
			imax = 6;
		} else {
			imax = 7;
		}
		// calculate and store most probable case
		logV[curr_state + 8 * t] = tomax[imax-6] +
										(est_emi[curr_state][(seq[t-1]-1) +
										emi_rows[curr_state] * (seq[t]-1)]);
		psi[curr_state][t] = imax;
		phi[curr_state][t] = 1;
	}

	free(seq);

	// TERMINATION
	t = seqLen; // i.e. an element after the sequence

	// initiate the termination to inf
	for (int j = 0; j < S - 1; j++) {
		logV[j + 8 * t] = -INFINITY;
	}

	// must terminate in null state, state of last sequence element must be
	// either R' or null
	curr_state = 7;
	tomax[0] = logV[6 + 8 * (t-1)] + (est_trans[6 * 8 + curr_state]);
	tomax[1] = logV[7 + 8 * (t-1)] + (est_trans[7 * 8 + curr_state]);
	if (tomax[0] >= tomax[1]) {
		imax = 6;
	} else {
		imax = 7;
	}

	// assign most probable case (R' or null in prev state)
	logV[curr_state + 8 * t] = tomax[imax-6];
	psi[curr_state][t] = imax;
	phi[curr_state][t] = 1;

	/*
	calculate the most probable path of the entire sequence by backtracking from
	the end to the begininning. The most probable path is the path with largest logV
	*/

	// begin path search from the terminating element
	int jmax = 0;
	path[t] = 7;
	d = phi[7][t];
	t = seqLen - d; // the last element of the previous state

	while (t >= 0) {
		jmax = psi[path[t+d]][t+d]; // most probable prev state

		if (t > 0) { // if backtracking not done
			// assign the corresponding state to the path
			d = phi[jmax][t];

			if (t - d + 1 < 0) { // the state covers the rest of the sequence
				for (int i = 0; i <= t; i++) {
					path[i] = psi[path[t+d]][t+d];
				}
				t = -1;
			} else {
				for (int i = t - d + 1; i <= t; i++) {
					path[i] = jmax;
				}
				t = t - d;
			}
		} else {
			path[t] = jmax;
			t = -1;
		}
	}

	for (int i = 0; i < S; i++) {
		free(psi[i]);
		free(phi[i]);
	}
	free(psi);
	free(phi);

}
