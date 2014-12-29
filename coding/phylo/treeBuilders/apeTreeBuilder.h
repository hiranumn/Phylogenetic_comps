/*
 * apeTreeBuilder.h
 * Makes a tree object according to the commonly accepted ape phylogeny.
 * The data was taken from Perelman, Polina, et al. "A molecular phylogeny of living primates." PLoS genetics 7.3 (2011): e1001342.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */
#include "tree.h"

#ifndef ATB
#define ATB

class ApeTreeBuilder
{
 public:
  /*
   *returns a commonly accepted ape phylogeny tree. 
   */
  Tree makeTree();
};

#endif
