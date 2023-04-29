cmake_minimum_required (VERSION 2.8)
project (OversamplingOscillators)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)
include_directories(${PORTED_PATH}/plugins/)
include_directories(${PORTED_PATH}/plugins/ChowDSP)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

add_library(SawOS MODULE OversamplingOscillators.cpp)
add_library(SinOscOS MODULE OversamplingOscillators.cpp)
add_library(TriOS MODULE OversamplingOscillators.cpp)
add_library(VarSawOS MODULE OversamplingOscillators.cpp)
add_library(SquareOS MODULE OversamplingOscillators.cpp)
add_library(SawBL MODULE OversamplingOscillators.cpp)
add_library(SquareBL MODULE OversamplingOscillators.cpp)
add_library(ImpulseBL MODULE OversamplingOscillators.cpp)
add_library(TriBL MODULE OversamplingOscillators.cpp)

add_library(SawH MODULE OversamplingOscillators.cpp)