EXE=ocean
FC=gfortran
FFLAGS=-cpp -ffree-form -fdefault-real-8 -g -c 
LDFLAGS=-g
SRC=mod_data.f90 mpdata_adiff.f90 ocean.f90 

OBJ=$(patsubst %.f90,%.o,$(SRC))

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $^ -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm *.o *.mod
