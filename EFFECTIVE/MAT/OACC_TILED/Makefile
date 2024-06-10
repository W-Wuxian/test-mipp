CC = pgcc
CFLAGS = -O3 -acc=gpu -Minfo=all -gpu=cc60

all: oacc_tiled_matmul

oacc_matmul: oacc_tiled_matmul.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f oacc_tiled_matmul