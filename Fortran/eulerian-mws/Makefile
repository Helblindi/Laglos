OPT = -O4 -march=native -mtune=native -Ofast
#OPT = -Wall -fimplicit-none -fbounds-check
lambda: eulerian_arbitrary_eos.o main.o
	gfortran $(OPT) main.o eulerian_arbitrary_eos.o -o lambda
main.o: main.f90
	gfortran -c $(OPT) main.f90
eulerian_arbitrary_eos.o: eulerian_arbitrary_eos.f90
	gfortran -c $(OPT) eulerian_arbitrary_eos.f90
clean:
	rm *.o *.mod lambda
