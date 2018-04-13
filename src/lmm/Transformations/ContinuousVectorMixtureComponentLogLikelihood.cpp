/*
 *  ContinuousVectorMixtureComponentLogLikelihood.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 17/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <DAG/DAG.h>
#include <lmm/Transformations/ContinuousVectorMixtureComponentLogLikelihood.h>
#include <DAG/RandomVariable.h>

using namespace gcat;

namespace gcat_lmm {
	
const string ContinuousVectorMixtureComponentLogLikelihoodParameterNames[0];

ContinuousVectorMixtureComponentLogLikelihood::ContinuousVectorMixtureComponentLogLikelihood(string rv_name, string distribution_name, string name, DAG* dag) : DAGcomponent(name,dag,"ContinuousVectorMixtureComponentLogLikelihood"), Transformation(ContinuousVectorMixtureComponentLogLikelihoodParameterNames,0), _distribution_name(distribution_name) {
	_rv = getDAG()->get_random_variable(rv_name);
}

ContinuousVectorMixtureComponentLogLikelihood::ContinuousVectorMixtureComponentLogLikelihood(const ContinuousVectorMixtureComponentLogLikelihood& x) : DAGcomponent(x), Transformation(x), _distribution_name(x._distribution_name), _cvm(x._cvm), _rv(x._rv), _L(x._L) {
}

string ContinuousVectorMixtureComponentLogLikelihood::validate() const {
	// For starters, do not allow this to parameterize anything else
	if(n_child_distributions()+n_child_transformations()>0) {
		string errTxt = "ContinuousVectorMixtureComponentLogLikelihood: object " + name() + " may not parameterize other objects";
		return errTxt;
	}
	// Secondly, finish identifying the distribution to which it relates (Distributions
	// are initialized after Transformations, so this could not be done in the constructor)
	_cvm = dynamic_cast<ContinuousVectorMixture*>(getDAG()->get_distribution(_distribution_name));
	if(!_cvm) {
		string errTxt = "ContinuousVectorMixtureComponentLogLikelihood: object " + name() + " cannot find ContinuousVectorMixture object " + _distribution_name;
		return errTxt;
	}
	const ContinuousVectorVariable* ct = dynamic_cast<const ContinuousVectorVariable*>(_rv);
	if(!ct) {
		string errTxt = "ContinuousVectorMixtureComponentLogLikelihood: " + _rv->name() + " is not of type ContinuousVectorVariable";
		return errTxt;
	}
	_L = _cvm->nmix();
	return "";
}

int ContinuousVectorMixtureComponentLogLikelihood::length() const {
	return _L;
}

// NB: Repeatedly calling get_double() gives marginal draws, whereas get_doubles() draws a complete path
// For that reason, the print() method is overwritten
double ContinuousVectorMixtureComponentLogLikelihood::get_double(const int i) const {
	if(i<0) error("ContinuousVectorMixtureComponentLogLikelihood::get_double(): index cannot be negative");
	vector<mydouble> x = _cvm->stored_likelihoods(_rv);
	if(i>=_L) error("ContinuousVectorMixtureComponentLogLikelihood::get_double(): index too large");
	return x[i].LOG();
}

vector<double> ContinuousVectorMixtureComponentLogLikelihood::get_doubles() const {
	vector<mydouble> x = _cvm->stored_likelihoods(_rv);
	vector<double> ret(x.size());
	for(int i=0;i<x.size();i++) ret[i] = x[i].LOG();
	return ret;
}

bool ContinuousVectorMixtureComponentLogLikelihood::has_changed(const int i) const {
	return true;
}

vector<bool> ContinuousVectorMixtureComponentLogLikelihood::has_changed() const {
	return vector<bool>(length(),true);
}

bool ContinuousVectorMixtureComponentLogLikelihood::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		default:
			error("ContinuousVectorMixtureComponentLogLikelihood::check_parameter_type(): parameter not found");
	}
	return false;
}

/*void ContinuousVectorMixtureComponentLogLikelihood::print(ostream& out, string sep) {
	int i;
	vector<double> x = get_doubles();
	for(i=0;i<length();i++) {
		if(i>0) out << sep;
		out << x[i];
	}
}*/
	
} // namespace gcat_lmm
