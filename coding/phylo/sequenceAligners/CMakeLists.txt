SET(sequenceAligners_SRCS 
   clustalW.cpp
   muscle.cpp
   centerStar.cpp
  ../treeBuilders/neighborJoining.cpp 
  ../treeBuilders/tree.cpp)

ADD_LIBRARY(sequenceAligners SHARED ${sequenceAligners_SRCS})
TARGET_LINK_LIBRARIES(sequenceAligners util)
