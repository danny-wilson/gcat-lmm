/*  Copyright 2014 Daniel Wilson.
 *
 *  Centre.cpp
 *  Part of the lmm library.
 *
 *  The lmm library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The lmm library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the lmm library. If not, see <http://www.gnu.org/licenses/>.
 */
#include <lmm/Transformations/Centre.h>

namespace gcat_lmm {
	
	const string CentreVectorTransformParameterNames[1] = {"x"};
	
	CentreVectorTransform::CentreVectorTransform(string name, DAG* dag) : DAGcomponent(name,dag,"CentreVectorTransform"), Transformation(CentreVectorTransformParameterNames,1), _x(0), _x_prev(0), _has_changed(0), _recalculate(true), _init(true) {
	}
	
	CentreVectorTransform::CentreVectorTransform(const CentreVectorTransform& x) : DAGcomponent(x), Transformation(x), _x(x._x), _x_prev(x._x_prev), _has_changed(x._has_changed), _recalculate(x._recalculate), _init(x._init) {
	}
	
	int CentreVectorTransform::length() const {
		return get_vector()->length();
	}
	
	double CentreVectorTransform::get_double(const int i) const {
		if(_recalculate) recalculate();
		return _x[i];
	}
	
	vector<double> CentreVectorTransform::get_doubles() const {
		if(_recalculate) recalculate();
		return _x;
	}
	
	bool CentreVectorTransform::has_changed(const int i) const {
		if(_recalculate) recalculate();
		return _has_changed[i];
	}
	
	vector<bool> CentreVectorTransform::has_changed() const {
		if(_recalculate) recalculate();
		return _has_changed;
	}
	
	bool CentreVectorTransform::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0:	// vector
				return(dynamic_cast<ContinuousVectorVariable*>(parameter));
			default:
				error("CentreVectorTransform::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	void CentreVectorTransform::set_vector(ContinuousVectorVariable* vec) {
		set_parameter(0,(Variable*)vec);
	}
	
	ContinuousVectorVariable const* CentreVectorTransform::get_vector() const {
		return (ContinuousVectorVariable const*)get_parameter(0);
	}
	
	void CentreVectorTransform::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
		if(sgl==Variable::_ACCEPT) {
			_has_changed = vector<bool>(length(),false);
		}
		else if (sgl==Variable::_REVERT) {
			_x = _x_prev;
			_has_changed = vector<bool>(length(),false);
		}
		else if(sgl==Variable::_SET) {
			_recalculate = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_recalculate = true;
		}
		// Call default implementation, which is to call Variable::send_signal_to_children(sgl)
		Transformation::receive_signal_from_parent(v,sgl);
	}
	
	void CentreVectorTransform::recalculate() const {
		if(_recalculate) {
			int i;
			const int n = length();
			_x_prev = _x;
			_x = get_vector()->get_doubles();
			double mn = 0.0;
			for(i=0;i<n;i++) {
				mn += _x[i];
			}
			mn /= (double)n;
			for(i=0;i<n;i++) {
				_x[i] -= mn;
			}
			if(_x_prev.size()==n) {
				for(i=0;i<n;i++) {
					_has_changed[i] = _x[i]!=_x_prev[i];
				}
			} else {
				_has_changed = vector<bool>(n,true);
			}
			_recalculate = false;
		}
	}
	
}	// namespace gcat_lmm
