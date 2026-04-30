CXX      = g++
NVCC     = nvcc

# Imath flag potrebbe essere necessario includerlo nei path
# Se Imath e' installato globalmente, -lImath basta.
CXXFLAGS = -std=c++17 -Wall -O2 -fopenmp
# Aggiunto -arch=native per assicurare il supporto a __half sulla tua GPU
NVCCFLAGS= -std=c++17 -O2 -Xcompiler -fopenmp -arch=native 
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

TARGET_CPU = $(BUILD)/test_parser
TARGET_GPU = $(BUILD)/test_parser_gpu

all: $(TARGET_CPU) $(TARGET_GPU)

cpu: $(TARGET_CPU)

gpu: $(TARGET_GPU)

$(TARGET_CPU): $(SRCS_CPU)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS_CPU) -o $(TARGET_CPU) -lImath
	@echo "CPU build completata!"

$(TARGET_GPU): $(SRCS_GPU)
	@mkdir -p $(BUILD)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(SRCS_GPU) -o $(TARGET_GPU) -lImath
	@echo "GPU build completata!"

clean:
	rm -rf $(BUILD)/*
	@echo "Pulizia completata."