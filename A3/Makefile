############################
# Makefile for assignment 3
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = quicksort quicksort_copy

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

