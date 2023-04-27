############################
# Makefile for lab 2
# Author: Malin Kallen 2019
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = exercise1a exercise1b exercise2 cart_shift_NSEW

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

