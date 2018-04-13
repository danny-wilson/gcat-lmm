/*  Copyright 2014 Daniel Wilson.
 *
 *  RescaleSymmetricMatrix.cpp
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
#include <lmm/Transformations/RescaleSymmetricMatrix.h>

namespace gcat_lmm {
	
	const string RescaleSymmetricMatrixTransformParameterNames[2] = {"scalar","matrix"};
	
	RescaleSymmetricMatrixTransform::RescaleSymmetricMatrixTransform(string name, DAG* dag) : DAGcomponent(name,dag,"RescaleSymmetricMatrixTransform"), Transformation(RescaleSymmetricMatrixTransformParameterNames,2), _x(), _x_prev(), _any_has_changed(true), _has_changed(), _recalculate(true), _init(true) {
	}
	
	RescaleSymmetricMatrixTransform::RescaleSymmetricMatrixTransform(const RescaleSymmetricMatrixTransform& x) : DAGcomponent(x), Transformation(x), _x(x._x), _x_prev(x._x_prev), _any_has_changed(x._any_has_changed), _has_changed(x._has_changed), _recalculate(x._recalculate), _init(x._init) {
	}
	
	int RescaleSymmetricMatrixTransform::length() const {
		return get_inmatrix()->length();
	}
	
	double RescaleSymmetricMatrixTransform::get_double(const int i, const int j) const {
		if(_recalculate) recalculate();
		return _x.safe(i,j);
	}
	
	vector<double> RescaleSymmetricMatrixTransform::get_lower_triangle() const {
		if(_recalculate) recalculate();
		const double sca = get_scalar()->get_double();
		vector<double> ret = get_inmatrix()->get_lower_triangle();
		int i;
		for(i=0;i<ret.size();i++) {
			ret[i] *= sca;
		}
		return ret;
	}
	
	HalfMatrix<double> RescaleSymmetricMatrixTransform::get_matrix() const {
		if(_recalculate) recalculate();
		return _x;
	}
	
	bool RescaleSymmetricMatrixTransform::has_changed(const int i, const int j) const {
		if(_recalculate) recalculate();
		return _has_changed.safe(i,j);
	}
	
	bool RescaleSymmetricMatrixTransform::has_changed() const {
		if(_recalculate) recalculate();
		return _any_has_changed;
	}
	
	bool RescaleSymmetricMatrixTransform::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0:	// scalar
				return(dynamic_cast<ContinuousVariable*>(parameter));
			case 1:	// vector
				return(dynamic_cast<SymmetricMatrixVariable*>(parameter));
			default:
				error("RescaleSymmetricMatrixTransform::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	void RescaleSymmetricMatrixTransform::set_scalar(ContinuousVariable* sca) {
		set_parameter(0,(Variable*)sca);
	}
	
	void RescaleSymmetricMatrixTransform::set_inmatrix(SymmetricMatrixVariable* vec) {
		set_parameter(1,(Variable*)vec);
	}
	
	ContinuousVariable const* RescaleSymmetricMatrixTransform::get_scalar() const {
		return (ContinuousVariable const*)get_parameter(0);
	}
	
	SymmetricMatrixVariable const* RescaleSymmetricMatrixTransform::get_inmatrix() const {
		return (SymmetricMatrixVariable const*)get_parameter(1);
	}
	
	void RescaleSymmetricMatrixTransform::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
		if(sgl==Variable::_ACCEPT) {
			_has_changed = HalfMatrix<bool>(length(),false);
			_any_has_changed = false;
		}
		else if (sgl==Variable::_REVERT) {
			_x = _x_prev;
			_has_changed = HalfMatrix<bool>(length(),false);
			_any_has_changed = false;
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
	
	void RescaleSymmetricMatrixTransform::recalculate() const {
		int i,j;
		const double sca = get_scalar()->get_double();
		const int n = length();
		_x_prev = _x;
		_x = get_inmatrix()->get_matrix();
		for(i=0;i<n;i++) {
			for(j=0;j<=i;j++) {
				_x[i][j] *= sca;
			}
		}
		if(_x_prev.n()==n) {
			_any_has_changed = false;
			for(i=0;i<n;i++) {
				for(j=0;j<=i;j++) {
					_has_changed[i][j] = _x[i][j]!=_x_prev[i][j];
					if(_has_changed[i][j]) _any_has_changed = true;
				}
			}
		} else {
			_any_has_changed = true;
			_has_changed = HalfMatrix<bool>(n,true);
		}
		_recalculate = false;
	}
	
}	// namespace gcat_lmm
