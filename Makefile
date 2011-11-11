.SUFFIXES: .f90

FC := gfortran

FFLAGS := -O3 -march=native -g
LIBS := -lrt -lblas

matmul_bench: matmul_bench.o
	$(FC) -o matmul_bench matmul_bench.o $(LIBS)

matmul_bench.o:
	$(FC) $(FFLAGS) -c $*.f90
