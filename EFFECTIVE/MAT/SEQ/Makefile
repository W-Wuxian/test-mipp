CC = icx
CFLAGS = -Wall -O3 -march=native

all: seq_matmul

seq_matmul: seq_matmul.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f seq_matmul