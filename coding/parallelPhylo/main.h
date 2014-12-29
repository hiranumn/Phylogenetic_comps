//The header file for main.cpp
//phylo
#include <iostream>
#include "util/species.h"
#include "util/options.h"
#include "util/sequenceData.h"
#include "util/distanceMatrix.h"
#include "treeBuilders/maximumParsimonyParallel.h"
#include "treeBuilders/maximumLikelihoodClimb.h"
#include "treeBuilders/tree.h"

#include "sequenceAligners/abstractAligner.h"
#include "sequenceAligners/clustalWParallel.h"
#include "sequenceAligners/clustalWOMP.h"

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <boost/graph/graphviz.hpp>

#include <omp.h>
