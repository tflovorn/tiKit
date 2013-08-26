FC = gfortran
FFLAGS =

LDFLAGS = -L/usr/lib -llapack -lblas

all: main.out

eigen.mod eigen.o: eigen.f90
	$(FC) -c -o eigen.o eigen.f90 $(FFLAGS)

main.out: eigen.mod eigen.o main.f90
	$(FC) -o main.out main.f90 eigen.o $(FFLAGS) $(LDFLAGS) 
