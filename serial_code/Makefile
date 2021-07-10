#This is a makefile for oblique shock reflection simulation.
shock_wave: main.o init.o solve.o
	gcc -Wall -o shock_wave main.o init.o solve.o -lm
main.o: main.c init.h solve.h grid.h
	gcc -Wall -c main.c 
init.o: init.c init.h grid.h
	gcc -Wall -c init.c
solve.o: solve.c solve.h init.h grid.h
	gcc -Wall -c solve.c

.PHONY: clean
clean:
	rm *.o *.txt shock_wave