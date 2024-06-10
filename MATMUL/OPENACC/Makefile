CC = pgcc
CFLAGS = -Wall -O3 -acc -ta=tesla

all: main

openacc_matrix_multiply: main.c
    $(CC) $(CFLAGS) -o $@ $<

clean:
    rm -f main