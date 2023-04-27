############################
# Makefile for lab 4
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = search

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

