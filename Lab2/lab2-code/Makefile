############################
# Makefile for lab 2
# Author: Malin Kallen 2019
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = alltoall arrays codetodebug integral2d IO_gather ring wave-parallel wave-persistent

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

