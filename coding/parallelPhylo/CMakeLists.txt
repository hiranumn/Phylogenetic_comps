CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(PhylogeneticReconstruction)

#######################
#Uncomment these for running on lab machines!!!
SET(BOOST_ROOT /usr/local/boost)
#######################

FIND_PACKAGE(Boost 1.4 COMPONENTS graph thread)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


ADD_SUBDIRECTORY(util)
ADD_SUBDIRECTORY(treeBuilders)
ADD_SUBDIRECTORY(sequenceAligners)
ADD_EXECUTABLE(phylo main.cpp)

TARGET_LINK_LIBRARIES(phylo util sequenceAligners treeBuilders ${Boost_LIBRARIES})
