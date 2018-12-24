# Set up the flags you want in FFLAGS and compiler in FC

# flags for GNU 
FFLAGS = -O0
FC = gfortran

# Suffixes rules to control how .o and .mod files are treated
.SUFFIXES:
.SUFFIXES: .f90 .o .mod 
.f90.o:
	$(FC) $(FFLAGS)  -c $< 
.f90.mod:
	$(FC) $(FFLAGS)  -c $< 

# Objects (OBJ) to be used

OBJ = SerraLINE.f90 parms.o io_mod.o functions_mod.o

#Compiles SerraLINE

SerraLINE:	$(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@

parms.o parms.mod: parms.f90
	$(FC) $(FFLAGS) -c parms.f90

#Instruction to clean mods
clean:
	rm *.o *.mod SerraLINE

