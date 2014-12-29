SET(treeBuilders_SRCS
  tree.cpp 
  maximumLikelihoodClimb.cpp
  neighborJoining.cpp
  maximumParsimonyParallel.cpp)

ADD_LIBRARY(treeBuilders SHARED ${treeBuilders_SRCS})

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

TARGET_LINK_LIBRARIES(treeBuilders util)
