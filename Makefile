
# These are for the ia64 cluster
FC= g77
CC= gcc
FLINK= g77
CLINK= gcc
CFLAGS=
# FFLAGS = -O2 -w -mcmodel=large -Vaxlib
LIB_DIR= ${HOME}/lib



CLIB= -lm
LINK= ar
LFLAGS= -rlv

OBJNR= nrutil.c

clean:	
	rm -f *.o

SRC1=	gasdev.for src3ft.f snrsq.f scopy.f sscal.f
OBJ1=	gasdev.o   src3ft.o snrsq.o scopy.o sscal.o
futil:	$(SRC1)
	$(FC) $(FFLAGS) -c $(SRC1) && $(LINK) $(LFLAGS) \
	libfutil.a $(OBJ1)
	mv libfutil.a $(LIB_DIR)

SRC0=	ftwrite.c ftread.c vscal.c nrutil.c cwrite_.c i3tensor.c alloc.c file_control.c
OBJ0= 	ftwrite.o ftread.o vscal.o nrutil.o cwrite_.o i3tensor.o alloc.o file_control.o
cutil:	$(SRC0)
	$(CC) $(CFLAGS) -c $(SRC0) && $(LINK) $(LFLAGS) \
	libcutil.a $(OBJ0)
	mv libcutil.a $(LIB_DIR)

SRC2=	ftwrite.c ftread.c
OBJ2= 	ftwrite.o ftread.o
ftio:	$(SRC2)
	$(CC) $(CFLAGS) -c $(SRC2) && $(LINK) $(LFLAGS) \
	libftio.a $(OBJ2)
	mv libftio.a $(LIB_DIR)

