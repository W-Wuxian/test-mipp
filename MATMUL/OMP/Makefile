CC = gcc
CFLAGS = -Wall -O3 -march=native -fopenmp

all: main

main: main.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f main