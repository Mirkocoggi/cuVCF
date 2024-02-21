CXX = g++
CPPFLAGS = -O3 -std=c++17
LIBS = -lz

S2GA:
	mkdir -p bin/
	$(CXX) src/VCFparser_mt.cpp -o bin/VCFparser -fopenmp $(CPPFLAGS) $(LIBS)

all:
	$(S2GA)

clean:
	rm -rf bin/

.PHONY:all clean