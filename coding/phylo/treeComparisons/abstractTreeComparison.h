/* 
 * abstractTreeComparison.h
 * Contains the abstract tree comparison class. Tree comparisons are used
 * to judge how dissimilar trees are.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/18/2014
 */

#ifndef ABSTREECOMPARE
#define ABSTREECOMPARE

#include "../treeBuilders/tree.h"

class AbstractTreeComparison
{
public:

    /*
     * Given two trees, returns the distance between them. Wrapper function
     * that just calls getDistances in a special case.
     * args: const reference to tree one and tree 2
     * return: the distance between tree 1 and tree two.
     */
    virtual double getDistance(const Tree& T1, const Tree& T2)
    {
	Tree copy(T2);
	vector<Tree> otherVector;
	otherVector.push_back(copy);
	return this->getDistances(T1, otherVector)[0];
    }
    
    /*
     * Given a "root" tree and a vector of trees to compare against, return
     * a vector of doubles corresponding to the ordered pairwise distances between
     * the "root" tree and the given vector of trees
     * args: a "root" tree to compare against and a vector of trees to compare against the
     * root tree
     * returns: a vector of distances from each of the comparison trees to the root tree.
     * The length of the return vector is equal to the length of the input tree vector.
     */
    virtual vector<double> getDistances(const Tree& comparison, const vector<Tree> &trees) = 0;
};


#endif
