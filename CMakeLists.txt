cmake_minimum_required(VERSION 3.24)
project(cuVCF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CUVCF_CUDA AUTO CACHE STRING "Build the GPU backend: AUTO, ON or OFF")
set_property(CACHE CUVCF_CUDA PROPERTY STRINGS AUTO ON OFF)
option(CUVCF_SANITIZE "Build CPU targets with AddressSanitizer + UBSan" OFF)

find_package(ZLIB REQUIRED)
find_package(Imath CONFIG REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# ── CPU backend ──────────────────────────────────────────────────
set(CPU_DEPS ZLIB::ZLIB Imath::Imath OpenMP::OpenMP_CXX)

add_executable(VCFparser_cpu src/CPUVersion/VCFparser_mt_col.cpp)
target_link_libraries(VCFparser_cpu PRIVATE ${CPU_DEPS})

pybind11_add_module(CPUParser src/CPUVersion/Bindings.cpp src/CPUVersion/VCFparser_mt_col.cpp)
target_link_libraries(CPUParser PRIVATE ${CPU_DEPS})

if(CUVCF_SANITIZE)
  foreach(t VCFparser_cpu CPUParser)
    target_compile_options(${t} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(${t} PRIVATE -fsanitize=address,undefined)
  endforeach()
endif()

install(TARGETS CPUParser LIBRARY DESTINATION .)

# ── GPU backend (optional) ───────────────────────────────────────
set(CUVCF_WITH_CUDA OFF)
if(NOT CUVCF_CUDA STREQUAL "OFF")
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    set(CUVCF_WITH_CUDA ON)
  elseif(CUVCF_CUDA STREQUAL "ON")
    message(FATAL_ERROR "CUVCF_CUDA=ON but no CUDA compiler was found")
  else()
    message(STATUS "cuVCF: no CUDA compiler found, GPU backend skipped")
  endif()
endif()

if(CUVCF_WITH_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    execute_process(
      COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
      RESULT_VARIABLE smi_rc OUTPUT_VARIABLE smi_out ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(smi_rc EQUAL 0 AND smi_out)
      set(CMAKE_CUDA_ARCHITECTURES native)
    else()
      set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 90)
    endif()
  endif()
  message(STATUS "cuVCF: CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")

  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  find_package(CUDAToolkit REQUIRED)

  set(GPU_DEPS ZLIB::ZLIB Imath::Imath OpenMP::OpenMP_CXX CUDA::cudart)
  set(GPU_FLAGS -Xcompiler=-fopenmp $<$<CONFIG:Debug>:-G -lineinfo>)

  # main.cu #includes Parser.cu, which #includes Kernels.cu: do not list them.
  add_executable(VCFparser_gpu src/GPUVersion/main.cu)
  target_compile_options(VCFparser_gpu PRIVATE ${GPU_FLAGS})
  target_link_libraries(VCFparser_gpu PRIVATE ${GPU_DEPS})

  # Bindings.cpp includes CUDA headers; compile it with nvcc like the old Makefile.
  set_source_files_properties(src/GPUVersion/Bindings.cpp PROPERTIES LANGUAGE CUDA)
  pybind11_add_module(GPUParser NO_EXTRAS src/GPUVersion/Bindings.cpp src/GPUVersion/Parser.cu)
  target_compile_options(GPUParser PRIVATE ${GPU_FLAGS} -Xcompiler=-w)
  target_link_libraries(GPUParser PRIVATE ${GPU_DEPS})

  install(TARGETS GPUParser LIBRARY DESTINATION .)
endif()
