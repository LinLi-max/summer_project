# This is a makefile for oblique shock reflection simulation.
# options: yes, no
PROFILING = yes
CC = mpicc
CFLAGS = -Wall
LDFLAGS = -lmpi -lm -g
PREP = scorep
PREP_CFLAGS = -lpapi

OBJECTS = main.c init.c decomp1d.c solve.c

shock_wave: $(OBJECTS)
ifeq ($(PROFILING),no) 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ 
else
	$(PREP) $(CC) $(CFLAGS) $(LDFLAGS) $(PREP_CFLAGS) -o $@ $^
endif

.PHONY: clean test
test:
	mpirun -n 2 ./shock_wave
clean:
	rm *.o *.txt shock_wave
