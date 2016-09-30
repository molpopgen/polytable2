CXXFLAGS=-I. -I.. -std=c++11 -Wall -W

all: src/PolyTable.o test.o
	$(CXX) $(CXXFLAGS) -o test test.o src/PolyTable.o -lgsl -lgslcblas

test.o: Sequence/PolyTable.hpp

src/PolyTable.o:Sequence/PolyTable.hpp

clean:
	rm -f test *.o
