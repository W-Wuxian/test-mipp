CC = pgcc
CFLAGS = -O3 -acc=gpu -Minfo=all -gpu=cc60

all: oacc_matmul

oacc_matmul: oacc_matmul.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f oacc_matmul