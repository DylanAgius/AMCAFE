# #!/bin/sh

CPP = /opt/petsc/arch-linux2-c-opt/bin/mpic++
CPPFLAGS = -std=c++11 -DADIOS2_USE_MPI -isystem /home/kteferra/Documents/research/software/ADIOS2/include

#CPPOOPTTFLAGS = -O2 -Wall
CPPOOPTTFLAGS = -Wall -Werror
#CPPOOPTTFLAGS = -O3

CPPINCLUDE = -I /opt/petsc/arch-linux2-c-opt/include/ 
METISLIB = /opt/petsc/arch-linux2-c-opt/lib/libmetis.so
adioslib=-Wl,-rpath,/home/kteferra/Documents/research/software/ADIOS2/lib /home/kteferra/Documents/research/software/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /home/kteferra/Documents/research/software/ADIOS2/lib/libadios2_cxx11.so.2.6.0 -Wl,-rpath-link,/home/kteferra/Documents/research/software/ADIOS2/lib


#SOURCES = $(wildcard *.C)
SOURCES = main1.C Grid.C BasePlate.C TempField.C Partition.C SampleOrientation.C VoxelsCA.C 
SOURCES2 = main1Scale.C Grid.C BasePlate.C TempFieldScale.C Partition.C SampleOrientation.C VoxelsCA.C 
SOURCES3 = main1val1.C Grid.C BasePlateval1.C TempFieldval1.C Partition.C SampleOrientation.C VoxelsCAval1.C
OBJECTS = $(SOURCES:.C=.o)
OBJECTS2 = $(SOURCES2:.C=.o)
OBJECTS3 = $(SOURCES3:.C=.o)


cafe: $(OBJECTS)
	@rm -f $@
	$(CPP) -o $@ $^ $(METISLIB) $(adioslib)

cafeScale: $(OBJECTS2)
	@rm -f $@
	$(CPP) -o $@ $^ $(METISLIB) $(adioslib)

cafeval1: $(OBJECTS3)
	@rm -f $@
	$(CPP) -o $@ $^ $(METISLIB) $(adioslib)

clean:
	rm -f *.o

#--------------------------------------------------------
# pattern rules for creating objects
.SUFFIXES: .C  # define the suffixes

%.o : %.C
	$(CPP) -c $(CPPFLAGS) $(CPPOPTFLAGS) $(CPPINCLUDE) $< -o $@

# pattern rules for creating objects
#--------------------------------------------------------
