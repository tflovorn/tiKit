FC = gfortran
FFLAGS =

LDFLAGS =

all: main.out

eigen.mod eigen.o: eigen.f90
	$(FC) $(FFLAGS) -c -o eigen.o eigen.f90

main.out: eigen.mod eigen.o main.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -o main.out main.f90 eigen.o
