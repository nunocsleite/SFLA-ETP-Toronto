
# In order to use C++ 11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")


cmake_minimum_required(VERSION 2.8)

# For Release compilation
set(CMAKE_BUILD_TYPE Release)


# About this project
#-------------------------------
project(MOEA)


#message(${PROJECT_SRC_DIR}/include)

# Include directories
#-------------------------------
#include_directories(${PROJECT_SRC_DIR}/include)

#SET(BOOST_INCLUDE_DIR "D:/Projecto/PhD/Libraries/boost_1_52_0/")
#SET(PARADISEO__EO_INCLUDE_DIR "D:/Projecto/PhD/Libraries/ParadisEO-2.0/eo/src/")
#SET(PARADISEO__MO_INCLUDE_DIR "D:/Projecto/PhD/Libraries/ParadisEO-2.0/mo/src/")
#SET(PARADISEO__MOEO_INCLUDE_DIR "D:/Projecto/PhD/Libraries/ParadisEO-2.0/moeo/src/")

#SET(BOOST_INCLUDE_DIR ${PROJECT_SRC_DIR}/include/boost_1_52_0/)
#SET(PARADISEO__EO_INCLUDE_DIR ${PROJECT_SRC_DIR}/include/ParadisEO-2.0/eo/src/)
#SET(PARADISEO__MO_INCLUDE_DIR ${PROJECT_SRC_DIR}/include/ParadisEO-2.0/mo/src/)
#SET(PARADISEO__MOEO_INCLUDE_DIR ${PROJECT_SRC_DIR}/include/ParadisEO-2.0/moeo/src/)


#message(${BOOST_INCLUDE_DIR})
#message(${PARADISEO__EO_INCLUDE_DIR})
#message(${PARADISEO__MO_INCLUDE_DIR})
#message(${PARADISEO__MOEO_INCLUDE_DIR})

#include_directories(${BOOST_INCLUDE_DIR}; ${PARADISEO__EO_INCLUDE_DIR}; ${PARADISEO__MO_INCLUDE_DIR}; ${PARADISEO__MOEO_INCLUDE_DIR})


include_directories(${MOEA_SOURCE_DIR}/include/boost_1_52_0/)
include_directories(${MOEA_SOURCE_DIR}/include/ParadisEO-2.0/eo/src/)
include_directories(${MOEA_SOURCE_DIR}/include/ParadisEO-2.0/mo/src/)
include_directories(${MOEA_SOURCE_DIR}/include/ParadisEO-2.0/moeo/src/)

message(${MOEA_SOURCE_DIR})

#include_directories(external/lib)
include_directories(lib)
include_directories(exes)

# Libraries
#-------------------------------
#set(EO_LIBS "D:/Projecto/PhD/Libraries/ParadisEO-2.0-build/eo/lib")
#set(MOEO_LIBS "D:/Projecto/PhD/Libraries/ParadisEO-2.0-build/moeo/lib")

#set(EO_LIBS "/../../../../../Libraries/Linux/ParadisEO-2.0/eo/lib)
#set(MOEO_LIBS "../../../../../Libraries/Linux/ParadisEO-2.0/moeo/lib")

#set(EO_LIBS ${PROJECT_SRC_DIR}/external/lib/eo)
#set(MOEO_LIBS ${PROJECT_SRC_DIR}/external/lib/moeo)


set(REGEX_LIB ${MOEA_SOURCE_DIR}/external/lib/)
message(${REGEX_LIB})


#message(${EO_LIBS})
#message(${MOEO_LIBS})

#link_directories(${EO_LIBS}; ${MOEO_LIBS})
link_directories(${REGEX_LIB})


# Add subfolders
#-------------------------------
add_subdirectory(lib)
add_subdirectory(exes)
