SET(treeBuilders_SRCS
  tree.cpp 
  maximumLikelihoodClimb.cpp 
  maximumLikelihoodProg.cpp
  maximumParsimonyClimb.cpp
  maximumParsimonyProg.cpp
  neighborJoining.cpp
  monteCarlo.cpp
  apeTreeBuilder.cpp)

ADD_LIBRARY(treeBuilders SHARED ${treeBuilders_SRCS})

TARGET_LINK_LIBRARIES(treeBuilders util)
