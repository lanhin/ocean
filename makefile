EXE=ocean
FC=ifort
FFLAGS=-cpp -O3 -DINNERLOOP -qopenmp -c -qopt-report=3 -align array64byte -xMIC-AVX512
#-ffree-form -fdefault-real-8 -xCOMMON-AVX512 -xMIC-AVX512 -xCORE-AVX512
LDFLAGS=-O3 -qopenmp -xMIC-AVX512
SRC=mod_data.f90 mpdata_adiff.f90 ocean.f90

OBJ=$(patsubst %.f90,%.o,$(SRC))

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $^ -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm *.o *.mod *.optrpt $(EXE)
