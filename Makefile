CXX = g++
CPPFLAGS = -O3 -std=c++17
DEBUGFLAGS = -std=c++17
LIBS = -lz -lHalf

VARSTRUCT:
	mkdir -p bin/
	$(CXX) src/VCFparser_mt.cpp -o bin/VCFparser -fopenmp $(CPPFLAGS) $(LIBS)
#	$(CXX) src/VCFparser_mt.cpp -o bin/VCFparser -fopenmp $(CPPFLAGS) $(LIBS) -g
#poi lancialo con gcc --args
VARCOL:
	mkdir -p bin/
	$(CXX) src/VCFparser_mt_col.cpp -o bin/VCFparser -fopenmp $(CPPFLAGS) $(LIBS)
#	$(CXX) src/VCFparser_mt_col.cpp -o bin/VCFparser -fopenmp $(CPPFLAGS) $(LIBS) -g

#poi lancialo con gcc --args
DEBUG:
	$(CXX) src/VCFparser_mt_col.cpp -o bin/VCFparser -fopenmp $(DEBUGFLAGS) $(LIBS) -g
all:
	$(VARSTRUCT)
	$(VARCOL)

clean:
	rm -rf bin/

.PHONY:all cleanin