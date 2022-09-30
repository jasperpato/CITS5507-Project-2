IDIR = ../include
CC = gcc
CFLAGS = -I$(IDIR) -std=c11 -Wall -Werror -fopenmp # -Wno-unused-function -g

ODIR = obj
LDIR = ../lib

LIBS = 

_DEPS = bond.h cluster.h site.h stack.h util.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = bond.o cluster.o percolate.o site.o stack.o util.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: percolate

percolate: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm ./percolate