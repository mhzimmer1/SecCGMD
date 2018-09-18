FC = gfortran 
EXE = sec3D.tachus
FLAGS = -O3
MAIN = main.f90
SOURCES = sys_param.f90 rand_num.f90 iolib.f90 polymer.f90 forceField.f90 channel.f90 integrator.f90 analysis.f90
OBJ = ${SOURCES:.f90=.o}

#Suffix Rules
.SUFFIXES: .f90 .o
.f90.o :
	$(FC) -c -O3 $<

$(EXE): $(MAIN) $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(MAIN) -o $(EXE)

clean:
	rm -f *.out *.o *.mod
