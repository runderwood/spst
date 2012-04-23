srcdir=./src/
builddir=./build/

CC=gcc
CFLAGS=-std=gnu99 -Wall

test:
	$(CC) -I. -I$(srcdir) -lm -lgdal $(CFLAGS) $(srcdir)sp.c $(srcdir)spst.c -o $(builddir)spst

