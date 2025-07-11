# ────────────────────────────────────────────────────────────────
#  P Y T H O N   /   P Y B I N D 1 1
# ────────────────────────────────────────────────────────────────
PYINC        := $(shell python3-config --includes)          # -I.../python3.x
PYBIND11INC  := $(shell python3 -m pybind11 --includes)     # -I.../pybind11
PYLDFLAGS    := $(shell python3-config --ldflags)           # -lpython3.x …
PYEXT        := $(shell python3-config --extension-suffix)  # .cpython-*.so

# ────────────────────────────────────────────────────────────────
#  C O M P I L A T I O N   F L A G S
# ────────────────────────────────────────────────────────────────
OPENMP          = -fopenmp
CPPFLAGS        = -O3 -std=c++17
GPUFLAGS        = -O3 -std=c++17 -arch=sm_89
DEBUG_CXXFLAGS  = -O0 -g -std=c++17 -fsanitize=address,undefined
DEBUG_GPUFLAGS  = -G -g -O0 -lineinfo -std=c++17 -arch=sm_89

# Librerie comuni
LIBS = -lz -I/usr/include/Imath -lImath        # aggiungi -lHalf se serve

# Toolchain
CXX  = g++
NVCC = nvcc

# Cartella di output
BIN_DIR = bin
$(BIN_DIR):
	@mkdir -p $@

# ────────────────────────────────────────────────────────────────
#  E X E C U T A B L E S   ( C P U   S O L O )
# ────────────────────────────────────────────────────────────────
VARSTRUCT: | $(BIN_DIR)
	$(CXX) src/VCFparser_mt.cpp -o $(BIN_DIR)/VCFparser $(OPENMP) $(CPPFLAGS) $(LIBS)

VARCOL: | $(BIN_DIR)
	$(CXX) src/CPUVersion/VCFparser_mt_col.cpp -o $(BIN_DIR)/VCFparser $(OPENMP) $(CPPFLAGS) $(LIBS)

VARCOL_DEBUG: | $(BIN_DIR)
	$(CXX) src/CPUVersion/VCFparser_mt_col.cpp \
	       -o $(BIN_DIR)/VCFparser_debug \
	       $(OPENMP) $(DEBUG_CXXFLAGS) $(LIBS)

# ────────────────────────────────────────────────────────────────
#  E X E C U T A B L E S   ( G P U )
# ────────────────────────────────────────────────────────────────
GPU: | $(BIN_DIR)
	$(NVCC) $(GPUFLAGS) --ptxas-options=-v \
	       -o $(BIN_DIR)/VCFparser src/GPUVersion/main.cu -Xcompiler $(OPENMP)

DEBUG: | $(BIN_DIR)
	$(NVCC) $(DEBUG_GPUFLAGS) \
	       -o $(BIN_DIR)/VCFparser src/GPUVersion/main.cu -Xcompiler $(OPENMP)

# ────────────────────────────────────────────────────────────────
#  P Y T H O N   M O D U L E S
# ────────────────────────────────────────────────────────────────
# GPU bindings
PYBIND: | $(BIN_DIR)
	$(NVCC) $(GPUFLAGS) $(PYINC) $(PYBIND11INC) \
	       -Xcompiler "-fPIC -w $(OPENMP)" -shared \
	       -o $(BIN_DIR)/GPUParser.so \
	       src/GPUVersion/Bindings.cpp src/GPUVersion/Parser.cu \
	       $(LIBS) $(PYLDFLAGS)

# CPU bindings (release)
CPUBIND: | $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(PYINC) $(PYBIND11INC) -shared -fPIC \
	      -o $(BIN_DIR)/CPUParser.so \
	      src/CPUVersion/Bindings.cpp src/CPUVersion/VCFparser_mt_col.cpp \
	      $(OPENMP) $(LIBS) $(PYLDFLAGS)

# CPU bindings con AddressSanitizer + UBSan (debug)
CPUBIND_DEBUG: | $(BIN_DIR)
	$(CXX) $(DEBUG_CXXFLAGS) $(PYINC) $(PYBIND11INC) -shared -fPIC \
	      -o $(BIN_DIR)/CPUParser_debug$(PYEXT) \
	      src/CPUVersion/Bindings.cpp src/CPUVersion/VCFparser_mt_col.cpp \
	      $(OPENMP) $(LIBS) $(PYLDFLAGS)

# ────────────────────────────────────────────────────────────────
#  B U I L D   A L L   &   C L E A N
# ────────────────────────────────────────────────────────────────
all: VARSTRUCT VARCOL GPU PYBIND CPUBIND

clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean VARSTRUCT VARCOL GPU DEBUG PYBIND CPUBIND CPUBIND_DEBUG
