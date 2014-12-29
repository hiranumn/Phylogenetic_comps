SET(util_SRCS 
	      species.cpp
	      options.cpp
	      sequenceData.cpp
	      mutationMap.cpp
	      aminoAlphabet.cpp
	      aminoSub.cpp
	      distanceMatrix.cpp
	      )

ADD_LIBRARY(util SHARED ${util_SRCS}) 
