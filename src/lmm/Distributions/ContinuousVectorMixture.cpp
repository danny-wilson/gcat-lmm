/*  Copyright 2012 Daniel Wilson.
 *
 *  ContinuousVectorMixture.cpp
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
#include <lmm/Distributions/ContinuousVectorMixture.h>
#include <sstream>

using std::stringstream;

namespace gcat {
	
	const string ContinuousVectorMixtureParameterNames[1] = {"p"};

	const string* ContinuousVectorMixtureDistributionNames(const int n) {
		string* ret = new string[n];
		int i;
		for(i=0;i<n;i++) {
			stringstream s;
			s << "distribution" << i;
			ret[i] = s.str();
		}
		return ret;
	}
	
	ContinuousVectorMixture::ContinuousVectorMixture(const int length, const int n, string name, DAG* dag) : DAGcomponent(name,dag,"ContinuousVectorMixture"), CompoundDistribution(ContinuousVectorMixtureDistributionNames(n),n,ContinuousVectorMixtureParameterNames,1), _length(length), _n(n), _last_signal(_INITIALIZE) {
	}
	
	ContinuousVectorMixture::ContinuousVectorMixture(const ContinuousVectorMixture& x) : DAGcomponent((const DAGcomponent&)x), CompoundDistribution((const CompoundDistribution&)x), _length(x._length), _n(x._n), _last_signal(x._last_signal) {
	}
	
	bool ContinuousVectorMixture::check_random_variable_type(RandomVariable* random_variable) {
		ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(random_variable);
		// Check the length is compatible
		if(!(rv==0) && rv->length()!=length()) {
			stringstream errMsg;
			errMsg << "ContinuousVectorMixture::check_random_variable_type(): ContinuousVectorVariable object " << rv->name() << " has length " << rv->length();
			errMsg << " whereas ContinuousVectorMixture object " << name() << " has length " << length();
			error(errMsg.str().c_str());
		}
		return rv;
	}
	
	bool ContinuousVectorMixture::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0:	//	p
				return(dynamic_cast<ContinuousVectorVariable*>(parameter));
			default:
				error("ContinuousVectorMixture::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	void ContinuousVectorMixture::set_p(ContinuousVectorVariable* p) {
		set_parameter(0,(Variable*)p);
	}
	
	ContinuousVectorVariable const* ContinuousVectorMixture::get_p() const {
		return (ContinuousVectorVariable const*)get_parameter(0);
	}
	
	mydouble ContinuousVectorMixture::likelihood(const RandomVariable* rv, const Value* val) {
		_x = ((ContinuousVectorVariable*)val)->get_doubles();
		_p = get_p()->get_doubles();
		if(_p.size()!=_n) error("ContinuousVectorMixture::likelihood(): parameter p has wrong number of components");
		mydouble ret(0.0);
		_stored_likelihoods[rv] = vector<mydouble>(_n);
		for(int i=0;i<_n;i++) {
			_stored_likelihoods[rv][i] = get_parent(i)->likelihood(this,to_Value());
			ret += mydouble(_p[i]) * _stored_likelihoods[rv][i];
		}
		return ret;
	}
	
	// Get value at position i
	double ContinuousVectorMixture::get_double(const int i) const {
		if(!(i>=0 && i<_x.size())) error("ContinuousVectorMixture::get_double(): position out of range");
		return _x[i];
	}
	
	// Get vector of values
	vector<double> ContinuousVectorMixture::get_doubles() const {
		return _x;
	}
	
	int ContinuousVectorMixture::length() const {
		return _length;
	}
	
	// Indicate if variables have changed. Warning: inefficient
	bool ContinuousVectorMixture::has_changed(const int i) const {
		return true;
	}
	
	vector<bool> ContinuousVectorMixture::has_changed() const {
		return vector<bool>(_x.size(),true);
	}

	int ContinuousVectorMixture::nmix() const {
		return _n;
	}

	vector<mydouble> ContinuousVectorMixture::stored_likelihoods(const RandomVariable* rv) const {
		std::map< const RandomVariable*, vector<mydouble> >::const_iterator it = _stored_likelihoods.find(rv);
		if(it==_stored_likelihoods.end()) error("ContinuousVectorMixture::stored_likelihoods(): cannot find random variable");
		return it->second;
	}
	
	void ContinuousVectorMixture::receive_signal_from_parent(const Distribution* dist, const Signal sgl) {
		if(sgl!=_last_signal) {
			if(sgl==_SET) {
			}
			else if(sgl==_PROPOSE) {
				_previous_stored_likelihoods = _stored_likelihoods;
			}
			else if(sgl==_ACCEPT) {
			}
			else if(sgl==_REVERT) {
				_stored_likelihoods = _previous_stored_likelihoods;
			}
			else error("ContinuousVectorMixture::receive_signal_from_parent(): unexpected signal");
			_last_signal = sgl;
			propagate_signal_to_children(sgl);
		}
	}	
	
} // namespace gcat

