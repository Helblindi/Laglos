OPT = -O4 -march=native -mtune=native -Ofast
#OPT = -Wall -fimplicit-none -fbounds-check
lambda: vdw.o main_vdw.o
	gfortran $(OPT) main_vdw.o vdw.o -o lambda
main_vdw.o: main_vdw.f90
	gfortran -c $(OPT) main_vdw.f90
vdw.o: vdw.f90
	gfortran -c $(OPT) vdw.f90
clean:
	rm *.o *.mod lambda