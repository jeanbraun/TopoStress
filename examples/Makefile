.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .f95 .e .r .y .yr .ye .l .p .sh .csh .h

F90 = gfortran

.f90.o:
	$(F90) -c -O $*.f90

all: $(TOPOSTRESSES_OBJ)
	make RUN

RUN:
	$(F90) -O TopoStresses_run.f90 vtk.f90 -o TopoStresses_run.exe -L../lib -lTopoStresses
