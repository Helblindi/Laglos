OPT = -O4 -march=native -mtune=native -Ofast
#OPT = -Wall -fimplicit-none -fbounds-check
lambda: lagrangian_arbitrary_eos.o main.o
	gfortran $(OPT) main.o lagrangian_arbitrary_eos.o -o lambda
main.o: main.f90
	gfortran -c $(OPT) main.f90
lagrangian_arbitrary_eos.o: lagrangian_arbitrary_eos.f90
	gfortran -c $(OPT) lagrangian_arbitrary_eos.f90
.PHONY : clean
clean :
	rm *.o *.mod $(objects) lambda
