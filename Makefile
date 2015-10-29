FC = gfortran

PROGRAM = run.exe
EXE0 = p__ehy.o e__ehy.o
PRM0 = c__eprm.o v__in.o
OBJ0 = p__initmg.o p__init.o p__passv.o p__passx.o
RAND0 = p__mtrnd.o

all : $(PROGRAM)

$(PRM0): %.o : %.f90
	$(FC) -c $<
$(RAND0): %.o : %.f90
	$(FC) -c $<
$(OBJ0): %.o : %.f90
	$(FC) -c $<
$(EXE0): %.o : %.f90
	$(FC) -c $<

$(PROGRAM): $(PRM0) $(RAND0) $(OBJ0) $(EXE0)
	$(FC) -o $@ $(EXE0) $(PRM0) $(RAND0) $(OBJ0)

clean:
	rm -f *.o *.mod
