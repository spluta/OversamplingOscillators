cmake_minimum_required (VERSION 2.8)
project (OversamplingOscillators)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)
include_directories(/Users/spluta1/Downloads/portedplugins-main/plugins/)
include_directories(/Users/spluta1/Documents/SC/SuperColliderBuild/portedplugins/plugins/ChowDSP)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

add_library(OSaw2 MODULE OversamplingOscillators.cpp)
add_library(OVarSaw2 MODULE OversamplingOscillators.cpp)