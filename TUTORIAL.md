HattCI tutorial
====

HattCI is a C-program for the identification of attC sites in any type of DNA data. It uses a hidden Markov model (HMM) to describe each part of the attC site in a probabilistic manner.

### Download

HattCI 1.0b [4.6 Mb] can be downloaded [here](https://github.com/maribuon/HattCI/archive/master.zip "Download HattCI 1.0b") .

### Installation

To install HattCI, save HattCI-master.zip to *yourdir*, enter a terminal, and type:

 $ cd *yourdir* <br>
 $ unzip HattCI-master.zip <br>
 $ cd HattCI-master<br>
 $ make<br>

Make sure *yourdir* is on your path, so that hattci.out is available.

### Usage

For the simplest way to run HattCI, start a terminal and then: <br>

$ hattci.out fasta_file output_file

### Arguments

-b: <br>
Processes both the ordinary sequences and the complementary sequence. <br>
To run for both directions, type:

$ hattci.out -b fasta_file  output_file<br><br>

-s x: <br>
HattCI reads x sequences at a time and processes them before reading the next x sequences, in order to avoid overextending RAM. Default is 1000 sequences. This flag gives the option to manually choose number of sequences to read, in the case of large sequences. <br>
To specify the chunk size when reading sequences, type for instance: 

 $ hattci.out -s 100000  fasta_file  output_file <br><br>

-t x:<br>
HattCI may run a large part of the computations in parallel, i.e. let different threads process a set of sequences, which in turn gives a reduced computation time. The parallelization works best when processing larger chunks of sequences at a time. Default number of threads are 1.<br>
To specify the number of threads, type for instance

$ hattci.out -t 4  fasta_file  output_file<br><br>

### Examples

#### 1. Annotating integrons

This example is about how to use HattCI to annotate integrons. We will use a small subset of our reference dataset: 3 sequences containg an integron with 3 gene cassettes, and therefore 3 attC sites each. The fasta file can be downloaded [here [16 kb]](https://github.com/maribuon/HattCI/blob/master/data/tutorial.fasta "Download data for example 1.").

In this case, it is not necessary to run on both strands, since we will assume the integrons have been correctly assembled so that the search is performed on the top strand by default. The files are also quite small, so we do not need to use more than one thread. We can thus call HattCI as follow:

$ hattci.out  tutorial.fasta   tutorial.out

Mbases processed: 0.015136<br>
Sequences read: 3 <br>
Sequences processed: 3 <br>
Hits: 9 <br>

The results are saved to tutorial.out. The results are divided in two parts. The first part is a table as below, where the first column is just a hit number, the second column is the accession number, start and end are the coordinates for the putative attC site, Vscore is the Viterbi score for the corresponding hit, and the following columns are the start positions of each attC site part. <br>

| hit | ID       |start|end |Vscore|	R'' |	spacer''|L''	 |loop | L' |spacer'|R'  |   |
| --- | ---      |---  |--- |---:  | --- | ---     | --- |---  | ---| ---   |--- |---|
| 1   |AB113580.1|2203 |2312|21.1  |2203 |2210     |2215 |2223 |2293|2300   |2306| 0 |
| 2   |AB113580.1|3353 |3470|5.7   |3353 |3360     |3365 |3373 |3452|3459   |3464| 0 |
| 3   |AB113580.1|4047 |4124|4.1   |4047 |4054     |4059 |4067 |4105|4112   |4118| 0 |
| 4   |AY183453.1|2603 |2662|16.7  |2603 |2610     |2615 |2623 |2643|2650   |2656| 0 |
| 5   |AY183453.1|3975 |4031|6.3   |3975 |3982     |3987 |3995 |4013|4020   |4025| 0 |
| 6   |AY183453.1|4825 |4884|20.6  |4825 |4832     |4837 |4845 |4865|4872   |4878| 0 |
| 7   |EU886977.1|1105 |1213|7.8   |1105 |1112     |1117 |1125 |1195|1202   |1207| 0 |
| 8   |EU886977.1|1433 |1492|12.6  |1433 |1440     |1445 |1453 |1473|1480   |1486| 0 |
| 9   |EU886977.1|2289 |2348|22.3  |2289 |2296     |2301 |2309 |2329|2336   |2342| 0 |

The second part of the results file, contain, for each hit, the DNA sequence, the annotation, where R corresponds to R'' and R', s corresponds to nucleotides found in the spacers, L is for L'' and L' and the dots (..) correspond to the loop. Below, we have as example the output for the first hit:

gtctaacaattcgttcaagccgacgttgcttcgtggcggcgcttgcgtgctacgctaagcttcgcacgccgcttgccactgcgcaccgcggcttaactcaggtgttaggg<br>
RRRRRRRsssssLLLLLLLL......................................................................LLLLLLLssssssRRRRRRR<br>
ID: gi|40645552|dbj|AB113580.1|

start position: 2203<br>
end position: 2312<br>
diffV: 21.123385<br>
hit #: 1<br>

For comparison, the coordinates found at GenBank for the sequences in our tutorial.fasta can be found [here [4 kb]](https://github.com/maribuon/HattCI/blob/master/data/tutorial.out "Download coordinates for data in example 1.").

#### 2. Quantifying attC sites in metagenomes

HattCI can be used to analyze large datasets. As an example we can use a <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/61/CAM_PROJ_SargassoSea.read.fa" target="_blank">metagenome from the Sargasso Sea [265 Mb]</a>.<br>

In this case, we want HattCI to analyze both strands, since the fragments have not been assembled. For this we will use the flag -b. In addition, we want to run parallel threads, since the file is considerably large. For this, we will use the flag -t 4. Then HattCI can be called as:

$ hattci.out  -b  -t 4  CAM_PROJ_SargassoSea.read.fa  sargasso.out

Mbases processed: 126.560981<br>
Sequences read: 606285<br>
Sequences processed: 606285<br>
Hits (top strand): 3608<br>
Hits (bottom strand): 4401<br>

This resulted in a total of 8009 hits (3608 top strand + 4401 bottom strand). The complete output file is here [3.2 Mb]. If we select only hits with Viterbi score (Vscore) > 7.5, we can see that there are no hits in the Sargasso sea.
