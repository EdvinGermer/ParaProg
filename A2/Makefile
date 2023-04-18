############################
# Makefile for assignment 2
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = matmul

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

