#This is a makefile for oblique shock reflection simulation.

shock_wave: main.o init.o decomp1d.o solve.o
	mpicc -Wall -o shock_wave main.o init.o decomp1d.o solve.o -lm 
main.o: main.c init.h solve.h grid.h
	mpicc -Wall -c main.c 
init.o: init.c init.h grid.h
	mpicc -Wall -c init.c
decomp1d.o: decomp1d.c decomp1d.h grid.h
	mpicc -Wall -c decomp1d.c
solve.o: solve.c solve.h init.h grid.h
	mpicc -Wall -c solve.c

.PHONY: clean test
test:
	mpirun -n 4 ./shock_wave
clean:
	rm *.o *.txt shock_wave
