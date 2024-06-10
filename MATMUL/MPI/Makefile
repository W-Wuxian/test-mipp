CC = mpicc
CFLAGS = -Wall -O3 -march=native

all: main

main: main.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f main