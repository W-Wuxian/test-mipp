CC = gcc
CFLAGS = -Wall -O3 -march=native -fopenmp

all: omp_matmul

omp_matmul: omp_matmul.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f omp_matmul