/*
 *  ContinuousVectorMixtureComponentLogLikelihood.h
 *  gcat
 *
 *  Created by Daniel Wilson on 17/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _CONTINUOUS_VECTOR_MIXTURE_COMPONENT_LOGLIKELIHOOD_H_
#define _CONTINUOUS_VECTOR_MIXTURE_COMPONENT_LOGLIKELIHOOD_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <DAG/Transformation.h>
#include <lmm/Distributions/ContinuousVectorMixture.h>

using namespace gcat;

namespace gcat_lmm {
	
class ContinuousVectorMixtureComponentLogLikelihood : public ContinuousVectorVariable, public Transformation {
private:
	string _distribution_name;					// Initially store just the name, until validate() is called
	mutable ContinuousVectorMixture* _cvm;		// Distribution
	RandomVariable* _rv;						// RV index wrt the distribution
	mutable int _L;								// Length: this is the number of components in the ContinuousVectorMixture distribution

public:
	// Constructor
	ContinuousVectorMixtureComponentLogLikelihood(string rv_name, string distribution_name, string name="", DAG* dag=0);
	// Copy constructor
	ContinuousVectorMixtureComponentLogLikelihood(const ContinuousVectorMixtureComponentLogLikelihood& x);
	
	// Implementation of virtual functions inherited from base classes
	// Get number of mixture components
	int length() const;
	// Get value at position i
	double get_double(const int i) const;
	// Get vector of values
	vector<double> get_doubles() const;
	// Has the value changed at position i?
	bool has_changed(const int i) const;
	// Has the value changed at each position?
	vector<bool> has_changed() const;
	// Type-checking for parameter(s)
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Overload method inherited from Transformation
	//void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	//void recalculate() const;

	// Overload method inherited from ContinuousVectorVariable
	//void print(ostream& out, string sep);

protected:
	// Overload method inherited from Component
	string validate() const;
};
	
} // gcat_lmm

#endif // _CONTINUOUS_VECTOR_MIXTURE_COMPONENT_LOGLIKELIHOOD_H_


