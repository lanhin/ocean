EXE=ocean
FC=ifort
FFLAGS=-fpp -O3 -FR -r8 -c -DINNERLOOP -DLENGTHARRAY -qopenmp -qopt-report=3 -align array64byte -fp-model fast=2 -heap-arrays -xMIC-AVX512
#-ffree-form -fdefault-real-8 -xCOMMON-AVX512 -xMIC-AVX512 -xCORE-AVX512 -qopt-streaming-stores always
LDFLAGS=-O3 -qopenmp -heap-arrays -xMIC-AVX512
SRC=mod_data.f90 mpdata_adiff.f90 ocean.f90

OBJ=$(patsubst %.f90,%.o,$(SRC))

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $^ libverify.a -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm *.o *.mod *.optrpt $(EXE)
