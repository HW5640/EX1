ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) 
CXXFLAGS  += $(ROOTCFLAGS)
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3

ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
P5640FLAGS  = -L${P5640LIB}/lib -lP5640  -I${P5640LIB}
GSLFLAGS    = -I${EBROOTGSL}/include/gsl  -I/usr/include/gsl -lgsl -lgslcblas


all: run1.1 run1.1gsl run1.2

run1.1:  run1.1.cpp ${P5640LIB}/lib/libP5640.a
	$(GXX) -o run1.1 run1.1.cpp $(ROOTFLAGS) $(P5640FLAGS)

run1.1gsl: run1.1gsl.cpp
	$(GXX) -orun1.1gsl run1.1gsl.cpp $(GSLFLAGS) $(ROOTFLAGS)

run1.2:  run1.2.cpp ${P5640LIB}/lib/libP5640.a
	$(GXX) -o run1.2 run1.2.cpp $(ROOTFLAGS) $(P5640FLAGS)


clean:
	rm -f run1.1 run1.1gsl run1.2 *dat *root
