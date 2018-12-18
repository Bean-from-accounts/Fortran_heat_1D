FC = gfortran
OPT = -g -O0 -fbounds-check

OBJ = module.o main.o subroutines.o

chaleur :	$(OBJ)
	$(FC) $(OPT) $(OBJ) -o chaleur

module.o :	module.f90
	$(FC) $(OPT) module.f90 -c

main.o :	main.f90
	$(FC) $(OPT) main.f90 -c

subroutines.o :	subroutines.f90
	$(FC) $(OPT) subroutines.f90 -c

clean :
	/bin/rm -f $(OBJ) *.mod

