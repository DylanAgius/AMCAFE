# #!/bin/sh

CPP = /home/dylana/Software/openmpi_4.1.5/bin/mpic++

#CPP = /sw/software/OpenMPI/4.0.5-gcccuda-2020b/mpic++
#CPPFLAGS = -std=c++11 -DADIOS2_USE_MPI -isystem /home/kteferra/Documents/research/software/ADIOS2/include
CPPFLAGS = -DADIOS2_USE_MPI -isystem /home/dylana/Software/ADIOS2_built/include

CPPOPTFLAGS = -g -O0
#CPPOPTFLAGS = -Wall
#CPPOPTFLAGS = -O3

CPPINCLUDE = -I/home/dylana/Software/metis_5.1.0/include 
METISLIB = /home/dylana/Software/metis_5.1.0/lib/libmetis.so


#adioslib=-Wl,-rpath,/home/kteferra/Documents/research/software/ADIOS2/lib /home/kteferra/Documents/research/software/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /home/kteferra/Documents/research/software/ADIOS2/lib/libadios2_cxx11.so.2.6.0 -Wl,-rpath-link,/home/kteferra/Documents/research/software/ADIOS2/lib


#adioslib=-Wl,-rpath,/usr/local/ADIOS2/lib /usr/local/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_cxx11.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_core_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_core.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_evpath.so /usr/local/ADIOS2/lib/libadios2_ffs.so.1.6.0 /usr/local/ADIOS2/lib/libadios2_atl.so.2.2.1 /usr/local/ADIOS2/lib/libadios2_dill.so.2.4.1 /usr/local/ADIOS2/lib/libadios2_enet.so.1.3.14 /usr/local/ADIOS2/lib/libadios2_taustubs.so /usr/local/ADIOS2/lib/libadios2_c_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_c.so.2.6.0 -Wl,-rpath-link,/usr/local/ADIOS2/lib 


adioslib=-Wl,-rpath,/home/dylana/Software/ADIOS2_built/lib /home/dylana/Software/ADIOS2_built/lib/libadios2_cxx11_mpi.so.2.9.0 /home/dylana/Software/ADIOS2_built/lib/libadios2_cxx11.so.2.9.0 -Wl,-rpath-link,/home/dylana/Software/ADIOS2_built/lib


extralib=-L/home/dylana/Software/libfabric/lib/libfabric.so.1
extraextralib=-L/usr/lib/libucp.so
libuct=-L/usr/lib/libuct.so
libucs=-L/usr/lib/libucs.so
libucm=-L/usr/lib/libucm.so 


#SOURCES = $(wildcard *.C)
SOURCES = main1.C Grid.C BasePlate.C TempField.C Partition.C SampleOrientation.C VoxelsCA.C 
OBJECTS = $(SOURCES:.C=.o)


cafe: $(OBJECTS)
	@rm -f $@
	$(CPP) -o $@ $^ $(adioslib) $(METISLIB) $(extralib) $(extraextralib) $(libuct) $(libucs) $(libucm)

clean:
	rm -f *.o

#--------------------------------------------------------
# pattern rules for creating objects
.SUFFIXES: .C  # define the suffixes

%.o : %.C
	$(CPP) -c $(CPPFLAGS) $(CPPOPTFLAGS) $(CPPINCLUDE) $< -o $@

# pattern rules for creating objects
#--------------------------------------------------------
