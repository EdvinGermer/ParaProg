############################
# Makefile for project
############################

CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = montecarlo
OBJS = prop.o

montecarlo: montecarlo.c prop.o
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS)

prop.o: prop.c
	$(CC) $(CCFLAGS) -c -o $@ $<

clean:
	$(RM) $(BINS) $(OBJS)
