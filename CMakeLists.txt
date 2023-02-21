cmake_minimum_required(VERSION 3.23)
project(cuda_1 CUDA)

set(CMAKE_CUDA_STANDARD 14)

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries     (PkgConfig::FFTW)

SET(MYCOMPILE_FLAGS "-lfftw3 -lm")

add_executable(cuda_1 main.cu code.cu)

set_target_properties(cuda_1 PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(cuda_1 fftw3 m)

