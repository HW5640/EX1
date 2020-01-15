ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) 
CXXFLAGS  += $(ROOTCFLAGS)
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3

ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
P5640FLAGS  = -L${P5640LIB}/lib -lP5640  -I${P5640LIB}
GSLFLAGS    = -I${EBROOTGSL}/include/gsl  -I/usr/include/gsl -lgsl -lgslcblas


all: RKnTest RKnDemo gsl_test duffing

RKnTest:  RKnTest.cpp ${P5640LIB}/lib/libP5640.a
	$(GXX) -o RKnTest RKnTest.cpp $(ROOTFLAGS) $(P5640FLAGS)

RKnDemo:  RKnDemo.cpp ${P5640LIB}/lib/libP5640.a
	$(GXX) -o RKnDemo RKnDemo.cpp $(ROOTFLAGS) $(P5640FLAGS)

gsl_test: gsl_test.cpp
	$(GXX) -ogsl_test gsl_test.cpp $(GSLFLAGS)

duffing:  duffing.cpp ${P5640LIB}/lib/libP5640.a
	$(GXX) -o duffing duffing.cpp $(ROOTFLAGS) $(P5640FLAGS)

clean:
	rm -f RKnDemo RKnTest gsl_test Projectile.png *~
