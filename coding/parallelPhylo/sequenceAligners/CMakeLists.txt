SET(sequenceAligners_SRCS clustalWParallel.cpp clustalWOMP.cpp)

ADD_LIBRARY(sequenceAligners SHARED ${sequenceAligners_SRCS})
TARGET_LINK_LIBRARIES(sequenceAligners util)
