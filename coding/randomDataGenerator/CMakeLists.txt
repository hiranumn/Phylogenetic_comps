CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(PhylogeneticReconstruction)

#######################
#Uncomment these for running on lab machines!!!
set(BOOST_ROOT /usr/local/boost)
#######################

FIND_PACKAGE(Boost 1.4 COMPONENTS graph)

######################
#Uncomment these for running on lab machines!!!
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARIY_DIR})
######################

ADD_LIBRARY(LIB ../phylo/util/species.cpp ../phylo/util/options.cpp ../phylo/util/mutationMap.cpp ../phylo/treeBuilders/tree.cpp ../phylo/util/sequenceData.cpp randomDataGenerator.cpp)
ADD_EXECUTABLE(randomDataGenerator  main.cpp)

TARGET_LINK_LIBRARIES(randomDataGenerator LIB ${Boost_LIBRARIES})
