CC = gcc
CFLAGS = -Wall -std=c99
ODIR = obj
LIBS = -lm

_OBJ = GaussianElimination.o backSubstitution.o solveTallMatrix.o nullSpaceCalculate.o Transform_To_RREF.o r_exchange.o printMatrix.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $^ $(CFLAGS)

compile_run: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)