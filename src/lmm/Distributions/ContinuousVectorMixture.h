/*  Copyright 2012 Daniel Wilson.
 *
 *  ContinuousVectorMixture.h
 *  Part of the gcat-core library.
 *
 *  The gcat-core library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The gcat-core library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the gcat-core library. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _CONTINUOUS_VECTOR_MIXTURE_H_
#define _CONTINUOUS_VECTOR_MIXTURE_H_
#include <DAG/CompoundDistribution.h>
#include <Variables/ContinuousVector.h>
#include <map>

using std::map;

namespace gcat {

class ContinuousVectorMixture : public ContinuousVectorVariable, public CompoundDistribution {
protected:
	// Length
	int _length;
	// Number of mixture components
	int _n;
	// Value of the vector
	vector<double> _x, _p;
	// Last accepted log-likelihoods: one vector per random variable
	map< const RandomVariable*, vector<mydouble> > _stored_likelihoods, _previous_stored_likelihoods;
	// Indicate whether things have changed
	Signal _last_signal;
public:
	// Constructor
	ContinuousVectorMixture(const int length, const int n, string name="", DAG* dag=0);
	// Copy constructor
	ContinuousVectorMixture(const ContinuousVectorMixture& x);

	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	bool check_parameter_type(const int i, Variable* parameter);
	void set_p(ContinuousVectorVariable* p);
	ContinuousVectorVariable const* get_p() const;
	
	// Compute log-likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	// Get value at position i
	double get_double(const int i) const;
	// Get vector of values
	vector<double> get_doubles() const;
	// Return length
	int length() const;
	// Indicate if variables have changed
	bool has_changed(const int i) const;
	vector<bool> has_changed() const;

	// Functions for ContinuousVectorMixtureComponentLogLikelihood transformation
	// Number of mixture components
	int nmix() const;
	// Last computed likelihoods
	vector<mydouble> stored_likelihoods(const RandomVariable* rv) const;
	// Redefine this function inherited from RandomVariable via CompoundDistribution. Interpret the signal then propagate to child RVs
	virtual void receive_signal_from_parent(const Distribution* dist, const Signal sgl);
	
};
	
} // namespace gcat

#endif // _CONTINUOUS_VECTOR_MIXTURE_H_