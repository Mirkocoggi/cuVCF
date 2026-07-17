CXX      = g++
NVCC     = nvcc

# Imath flag potrebbe essere necessario includerlo nei path
# Se Imath e' installato globalmente, -lImath basta.
CXXFLAGS = -std=c++17 -Wall -O2 -fopenmp
# Aggiunto -arch=native per assicurare il supporto a __half sulla tua GPU
# NVCC_CCBIN: host compiler per nvcc. Su alcune macchine (gcc/nvcc troppo recenti)
# serve forzare una versione piu' vecchia, es: make gpu NVCC_CCBIN=g++-12
# Lascia vuoto (make gpu NVCC_CCBIN=) per usare il default di nvcc.
NVCC_CCBIN ?= g++-12
# -I compat/nvcc_gcc_fix: nvcc non riconosce alcune intrinsics AMX/AVX512-BF16 dichiarate
# dagli header di sistema (bug del front-end di nvcc, non del progetto): vedi compat/nvcc_gcc_fix/README.md
# Innocuo se la macchina non ha questo bug.
NVCCFLAGS= -std=c++17 -O2 -Xcompiler -fopenmp -arch=native $(if $(NVCC_CCBIN),-ccbin $(NVCC_CCBIN)) -I compat/nvcc_gcc_fix
# -lineinfo keeps symbol/line info for Nsight Compute source correlation without disabling opts
NVCCFLAGS_PROF = $(NVCCFLAGS) -lineinfo
INCLUDES = -Iinclude -Isrc/CPUVersion -Isrc/GPUVersion
BUILD    = build

SRCS_CPU = src/CPUVersion/main_cpu.cpp \
           src/CPUVersion/CSVParser.cpp \
           src/CPUVersion/VCFReconstructor.cpp \
           src/Utils.cpp

SRCS_GPU = src/GPUVersion/main_gpu.cu \
           src/GPUVersion/VCFReconstructorGPU.cu \
           src/CPUVersion/CSVParser.cpp \
           src/Utils.cpp

TARGET_CPU  = $(BUILD)/test_parser
TARGET_GPU  = $(BUILD)/test_parser_gpu
TARGET_PROF = $(BUILD)/test_parser_gpu_prof

all: $(TARGET_CPU) $(TARGET_GPU)


cpu: $(TARGET_CPU)

gpu: $(TARGET_GPU)

$(TARGET_CPU): $(SRCS_CPU)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS_CPU) -o $(TARGET_CPU) -lImath
	@echo "CPU build completata!"

$(TARGET_GPU): $(SRCS_GPU)
	@mkdir -p $(BUILD)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(SRCS_GPU) -o $(TARGET_GPU) -lImath -lnvToolsExt
	@echo "GPU build completata!"

profile: $(SRCS_GPU)
	@mkdir -p $(BUILD)
	$(NVCC) $(NVCCFLAGS_PROF) $(INCLUDES) $(SRCS_GPU) -o $(TARGET_PROF) -lImath -lnvToolsExt
	@echo "Profile build completata → $(TARGET_PROF)"

clean:
	rm -rf $(BUILD)/*
	@echo "Pulizia completata."