##### HattCI: #####
HattCI is C-implementation for <i attC\i> site identification in DNA sequences described in:

HattCI: Fast and Accurate attC site Identification Using Hidden Markov Models

Pereira Mariana Buongermino, Wallroth Mikael, Kristiansson Erik, and Axelson-Fisk Marina. Journal of Computational Biology. November 2016, 23(11): 891-902. doi:10.1089/cmb.2016.0024. 

##### INSTALL: #####
1. Enter Terminal.
2. cd to the directory.
3. Type 'make' in terminal.

##### RUN: #####
To run HattCI, type
$ ./hattci.out <fasta file> <location to put the output file>
or
$ ./hattci.out the.fasta ./output/output.txt
Note that if output.txt already exists, it is overwritten.
Otherwise, it is created. If no output file is given, the output is written to
a standard file.

##### FLAGS: #####
HattCI handles options.
-b:
Processes both the ordinary sequences and the complementary
sequence.
To run for both directions, type
$ ./hattci.out -b <fasta file> <location to put the output file>

-s x:
HattCI reads x sequences at a time and processes them before reading the next x sequences, in order to avoid overextending RAM. Default is 1000 sequences.
This flag gives the option to manually choose number of sequences to read, in
the case of large sequences.
To specify the chunk size when reading sequences, type for instance
$ ./hattci.out -s 100000 <fasta file> <location to put the output file>

-t x:
HattCI may run a large part of the computations in parallel, i.e. let different threads process a set of sequences, which in turn gives a reduced computation time. The parallelization works best when processing larger chunks of sequences at a time. Default number of threads are 1.
To specify the number of threads, type for instance
$ ./hattci.out -t 4 <fasta file> <location to put the output file>

The flags works together, type for instance
$ ./hattci.out -b -s 100 <fasta file> <location to put the output file>

##### SCRIPT: #####
To run several fasta files, a script file is supplied.
The script processes all files in the input directory and prints results to the output directory, with one output file for each submitted input file.
A summary (i.e. output to terminal) of the processed data is printed to summary.txt.

To run script, type
$ ./runScript.sh
