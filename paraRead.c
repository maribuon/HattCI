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

Reads data from each file in the parameter folder
(i.e. ./parameters/<filename>.txt)
and stores the data in corresponding vector/matrix, altering input pointers.

Parameters:
iniPtr 		(input/output)	double pointer, initially null
transPtr 	(input/output)	double pointer, initially null
fPtr 		(input/output)	double pointer, initially null
emiPtr		(input/output) 	double pointer, initially null
Dmin		(input)			double array, min length of each state
Dmax		(input)			double array, max length of each state
Dmax_cum	(input)			double array, cumulative max length of each state
S			(input)			int, number of states
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void paraRead(double **iniPtr, double **transPtr, double ***fPtr,
			double ***emiPtr, int Dmin[], int Dmax[], int Dmax_cum[], int S) {
	FILE *file;
	int elem;
	double *temp;
	double *est_ini, *est_trans, **est_emi, **f;

	// read from est_ini.txt
	if ((file = fopen(DIR "/parameters/est_ini.txt", "r")) == NULL) {
		printf("Problems opening file est_ini.txt\n");
		exit(1);
	}
	est_ini = malloc(S * sizeof(double));
	elem = 0;
	while (fscanf(file, "%lf", (est_ini+elem)) != EOF) {
		*(est_ini + elem) = log(*(est_ini + elem));
		elem++;
	}
	fclose(file);

	// read from est_trans.txt
	if ((file = fopen(DIR "/parameters/est_trans.txt", "r")) == NULL) {
		printf("Problems opening file est_trans.txt\n");
		exit(1);
	}
	est_trans = malloc(S * S * sizeof(double));
	elem = 0;
	while (fscanf(file, "%lf", (est_trans+elem)) != EOF) {
		*(est_trans + elem) = log(*(est_trans + elem));
		elem++;
	}
	fclose(file);

	// read from est_emi.txt
	if ((file = fopen(DIR "/parameters/est_emi.txt", "r")) == NULL) {
		printf("Problems opening file f.txt\n");
		exit(1);
	}

	int divisor[8] = {256, 16, 32, 64, 28, 16, 256, 16}; // same as sizes
	est_emi = malloc(S * sizeof(double *));
	for (int i = 0; i < S; i++) {
		est_emi[i] = malloc(divisor[i] * sizeof(double));
	}
	for (int i = 1; i < 8; i++) {
		divisor[i] += divisor[i-1];
	}
	temp = malloc(divisor[S-1] * sizeof(double));
	elem = 0;
	while (fscanf(file, "%lf", (temp+elem)) != EOF) {
		elem++;
	}
	fclose(file);

	for (int j = 0; j < divisor[0]; j++) {
		est_emi[0][j] = log(temp[j]);
	}

	for (int i = 1; i < S; i++) {
		for (int j = divisor[i-1]; j < divisor[i]; j++) {
			est_emi[i][j-divisor[i-1]] = log(temp[j]);
		}
	}
	free(temp); // deallocate temporary data

	// read from f.txt
	if ((file = fopen(DIR "/parameters/f.txt", "r")) == NULL) {
		printf("Problems opening file f.txt\n");
		exit(1);
	}
	temp = malloc(Dmax_cum[S-1] * sizeof(double));
	elem = 0;
	while (fscanf(file, "%lf", (temp+elem)) != EOF) {
		elem++;
	}
	fclose(file);
	f = malloc(S * sizeof(double *));
	f[0] = malloc(Dmax_cum[0] * sizeof(double));
	for (int j = 0; j < Dmax_cum[0]; j++) {
		f[0][j] = log(temp[j]);
	}
	for (int i = 1; i < S; i++) {
		f[i] = malloc((Dmax_cum[i] - Dmax_cum[i-1]) * sizeof(double));
		for (int j = Dmax_cum[i-1]; j < Dmax_cum[i]; j++) {
			f[i][j - Dmax_cum[i-1]] = log(temp[j]);
		}
	}
	free(temp); // deallocate the read data from f
	// assign input pointers to corresponding values
	*iniPtr = est_ini;
	*transPtr = est_trans;
	*fPtr = f;
	*emiPtr = est_emi;
}
