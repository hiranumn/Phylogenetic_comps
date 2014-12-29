SET(util_SRCS 
	      species.cpp
	      options.cpp
	      sequenceData.cpp
	      mutationMap.cpp
	      distanceMatrix.cpp
	      parallelUtils.cpp
	      parallelUtilsOMP.cpp)



ADD_LIBRARY(util SHARED ${util_SRCS})
