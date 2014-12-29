SET(treeComparisons_SRCS
  pairwisePath.cpp
  quartet.cpp)

ADD_LIBRARY(treeComparisons SHARED ${treeComparisons_SRCS})

TARGET_LINK_LIBRARIES(treeComparisons util treeBuilders)
