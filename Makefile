FC := gfortran

FFLAGS := -O3 -march=native -g -ffpe-trap=invalid,zero,overflow
LIBS := -lrt -lblas

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

matmul_bench: matmul_bench.o
	$(FC) -o matmul_bench $< $(LIBS)

