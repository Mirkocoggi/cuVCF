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
