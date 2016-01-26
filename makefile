
FC=gfortran -cpp -Wall -pedantic -std=f2008
OPT= -O3

#for correct functionality it is required to put correct address of fgsl library

#PREFIX=/home/zbona/lib/fgsl
PREFIX=/usr/local
LIBS=  -L $(PREFIX)/lib -lfgsl
INCLUDE = -I /usr/include/fgsl


#DEFS=-DRAMP

DEFS=-DOKHRIM # use okhrimovsky scattering angle

#DEBUG= -g -ffpe-trap=invalid,zero,overflow,underflow,inexact  -fbacktrace  -fbounds-check
#DEBUG= -g0  -fbacktrace  -fbounds-check
#DEBUG=   -pg

reid.bin: main.o params.o ziggurat.o
	$(FC) $(OPT) $(DEFS)  $(DEBUG) $^ -o $@ $(LIBS)


main.o: main.f90 ziggurat.o params.o
	$(FC) $(OPT) $(DEFS) $(DEBUG)  -c $< -o $@ $(INCLUDE) $(LIBS)


params.o: params.f90
	$(FC) $(OPT) $(DEFS) $(DEBUG)  -c $< -o $@


ziggurat.o: ziggurat.f90
	$(FC) $(OPT) $(DEFS)  -c $< -o $@


clean:
	rm -fv *.mod *.o *.bin
