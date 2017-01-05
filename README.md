HattCI
====
![alt text](https://github.com/maribuon/HattCI/blob/master/hattci-logo.png "Logo")

HattCI is C-implementation for identification of *attC* sites in DNA sequences.

### Citation #####
HattCI: Fast and Accurate attC site Identification Using Hidden Markov Models <br>
Pereira, Mariana Buongermino; Wallroth, Mikael; Kristiansson, Erik; and Axelson-Fisk; Marina. Journal of Computational Biology. November 2016, 23(11): 891-902. [doi:10.1089/cmb.2016.0024](http://online.liebertpub.com/doi/abs/10.1089/cmb.2016.0024).

### License #####
HattCI is freely distributed under the GNU General Public License [(GPLv3)](https://opensource.org/licenses/GPL-3.0 "GNU General Public License version 3").

### Installation #####
1. Download HattCI, move it to a *yourDirectory* where you want it installed. Uncompress it.
2. Enter Terminal.
3. cd *yourDirectory*
4. Type 'make' in terminal.

Make sure *yourDirectory* is in your $PATH. 

### Usage #####
To run HattCI, type<br>
$ hattci.out \<fasta file> \<location to put the output file>

or <br>
$ hattci.out the.fasta ./output/output.txt

Note that if output.txt already exists, it is overwritten.<br>
Otherwise, it is created. If no output file is given, the output is written to
a standard file.

#### Options #####
-b:<br>
Processes both the ordinary sequences and the complementary
sequence.<br>
To run for both directions, type:<br>
$ hattci.out -b \<fasta file> \<location to put the output file>

-s x:<br>
HattCI reads x sequences at a time and processes them before reading the next x sequences, in order to avoid overextending RAM. Default is 1000 sequences.<br>
This flag gives the option to manually choose number of sequences to read, in the case of large sequences.<br>
To specify the chunk size when reading sequences, type for instance: <br>
$ hattci.out -s 100000 \<fasta file> \<location to put the output file>

-t x:<br>
HattCI may run a large part of the computations in parallel, i.e. let different threads process a set of sequences, which in turn gives a reduced computation time. The parallelization works best when processing larger chunks of sequences at a time. Default number of threads are 1.<br>
To specify the number of threads, type for instance:<br>
$ hattci.out -t 4 \<fasta file> \<location to put the output file>

The flags works together, type for instance: <br>
$ hattci.out -b -s 100 <fasta file> <location to put the output file>

### Uninstall #####
1. Open a Terminal
2. cd the directory where it is installed
3. make clean
4. delete the directory

