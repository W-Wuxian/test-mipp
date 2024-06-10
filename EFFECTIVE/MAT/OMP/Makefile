CC = icx
CFLAGS = -Wall -O3 -march=native -qopenmp

all: omp_matmul

omp_matmul: omp_matmul.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f omp_matmul