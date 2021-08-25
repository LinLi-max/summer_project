# Compilers and commands
CXX=		gcc
NVCC=		/usr/local/cuda/bin/nvcc
LINK=		/usr/local/cuda/bin/nvcc
DEL_FILE= 	rm -f

#Flags
CXXFLAGS	= -W -Wall -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops
NVCCFLAGS   = -O3 --use_fast_math --expt-relaxed-constexpr -arch sm_60 -lcudart

LIBS		= $(PARALLEL)
INCPATH		= /usr/local/cuda/include/

####### Files
SOURCES	= main.cu kernel.cu
OBJECTS	= main.o kernel.o

TARGET = shock_wave

all: main

main:$(OBJECTS)
	$(NVCC) $(OBJECTS) -o $(TARGET) $(NVCCFLAGS) -I$(INCPATH)

main.o: main.cu kernel.h
	$(NVCC) -c $(SOURCES) $(NVCCFLAGS) -I$(INCPATH) 

kernel.o: kernel.cu kernel.h
	$(NVCC) -c $(SOURCES) $(NVCCFLAGS) -I$(INCPATH)

.PHONY: clean
clean:
	$(DEL_FILE) $(OBJECTS) $(TARGET)
