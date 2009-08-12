/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Normalizer.h,v 1.13.2.1 2008/08/11 18:31:12 jasonw Exp $
 *******************************************************************************/
#ifndef NORMALIZER_H_
#define NORMALIZER_H_
class DataSet;
#include <set>
using namespace std;

namespace qranker {

class Normalizer
{
public:
	virtual ~Normalizer();
    virtual void setSet(set<DataSet *> & setVec, size_t numFeatures, size_t numRetentionFeatures){;}
    void normalizeSet(set<DataSet *> & setVec);
    void normalize(const double * in, double * out, size_t offset, size_t numFeatures);
    double normalize(const double in, size_t index) { return (in-sub[index])/div[index]; }
    virtual void unnormalizeweight(const vector<double>& in,vector<double>& out){;}
    virtual void normalizeweight(const vector<double>& in,vector<double>& out){;}
    static Normalizer * getNormalizer();
    static void setType(int type);
	const static int UNI = 0;
	const static int STDV = 1;
protected:
    Normalizer();
    static Normalizer * theNormalizer;
	static int subclass_type;
    size_t numFeatures, numRetentionFeatures;
    vector<double> sub;
    vector<double> div;
};

} // qranker namspace

#endif /*NORMALIZER_H_*/
