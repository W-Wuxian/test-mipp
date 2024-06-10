CC = icx
CFLAGS = -Wall -O3 -march=native

all: seq_tiled_matmul

seq_tiled_matmul: seq_tiled_matmul.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f seq_tiled_matmul