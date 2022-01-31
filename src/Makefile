.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .f95 .e .r .y .yr .ye .l .p .sh .csh .h

F90 = gfortran

TOPOSTRESSES_OBJ = TopoStresses_ctx.o TopoStresses_api.o TopoStress.o solve_skyline_renumber_n.o sloan.o

.f90.o:
	$(F90) -c -O $*.f90

all: $(TOPOSTRESSES_OBJ) LIB
	make LIB

LIB: $(FASTSCAPE_OBJ)
	ar rc ../lib/libTopoStresses.a $(TOPOSTRESSES_OBJ)
	f2py -c --f90flags=-O3 -m TopoStresses *.f90
	cp TopoStresses.cpython-*.so ../lib/TopoStresses.so
	cp TopoStresses.cpython-*.so ../examples/TopoStresses.so