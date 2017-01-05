OBJS = main.o paraRead.o hmm.o viterbi.o readFasta.o calcResults.o
CC = gcc
DIR = $(PWD)
CFLAGS = -O3 -fopenmp -g -std=c99 -D DIR='"$(PWD)"'

hattci.out: $(OBJS)
	$(CC)  -o $@ $(CFLAGS) $(OBJS) -lm

clean:
	rm -f $(OBJS) hattci.out
