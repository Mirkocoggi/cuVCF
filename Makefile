# Compilatore C++ e NVCC
CXX = g++
NVCC = nvcc

# Flags di compilazione per la versione ottimizzata (CPU e GPU)
CPPFLAGS = -O3 -std=c++17
GPUFLAGS = -O3 -std=c++17 -arch=sm_89
OPENMP = -fopenmp

# Flags di debug per la versione CPU (se necessario)
DEBUGFLAGS = -std=c++17
# Flag di debug per la versione GPU (cuda-gdb richiede -G e -g)
DEBUG_GPUFLAGS = -G -g -std=c++17 -arch=sm_89

# Librerie da linkare
LIBS = -lz -I/usr/include/Imath -lImath

# Target per compilare la versione CPU (VARSTRUCT)
VARSTRUCT:
	mkdir -p bin/
	$(CXX) src/VCFparser_mt.cpp -o bin/VCFparser $(OPENMP) $(CPPFLAGS) $(LIBS)

# Target per compilare la versione CPU (VARCOL)
VARCOL:
	mkdir -p bin/
	$(CXX) src/VCFparser_mt_col.cpp -o bin/VCFparser $(OPENMP) $(CPPFLAGS) $(LIBS)

# Target per compilare la versione GPU ottimizzata
GPU:
	mkdir -p bin/
	$(NVCC) $(GPUFLAGS) --ptxas-options=-v -o bin/VCFparser src/CuFiles/main.cu -Xcompiler $(OPENMP)

# Target DEBUG: compila la versione GPU per cuda-gdb con flag di debug (-G -g) runna con -> cuda-gdb --args ./bin/VCFparser -v data/IRBT3M.vcf -t 1

DEBUG:
	mkdir -p bin/
	$(NVCC) $(DEBUG_GPUFLAGS) -o bin/VCFparser src/CuFiles/main.cu -Xcompiler $(OPENMP)

# Target "all" compila tutti i target desiderati
all: VARSTRUCT VARCOL GPU

# Pulizia dei file compilati
clean:
	rm -rf bin/

.PHONY: all clean VARSTRUCT VARCOL GPU DEBUG
