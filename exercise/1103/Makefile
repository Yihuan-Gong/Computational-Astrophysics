#
#

F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

######

SOURCES=solver.f90 shooting1.f90 
 
CLEANSTUFF=rm -rf *.o *.mod shooting1 *~

all: solver shooting1 

solver: solver.f90
	$(F90) $(F90FLAGS) -c $^

shooting1: solver.f90 shooting1.f90
	$(F90) $(F90FLAGS) solver.f90 shooting1.f90 -o shooting1

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
