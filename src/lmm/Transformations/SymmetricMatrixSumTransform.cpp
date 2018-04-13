/*  Copyright 2014 Daniel Wilson.
 *
 *  SymmetricMatrixSumTransform.cpp
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
#include <lmm/Transformations/SymmetricMatrixSumTransform.h>
#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::stringstream;
using std::vector;

namespace gcat_lmm {
	
	const string* SymmetricMatrixSumTransformParameterNames(const int n) {
		string* ret = new string[n];
		int i;
		for(i=0;i<n;i++) {
			stringstream s;
			s << "operand" << i;
			ret[i] = s.str();
		}
		return ret;
	}
	
	SymmetricMatrixSumTransform::SymmetricMatrixSumTransform(const int n, string name, DAG* dag) : DAGcomponent(name,dag,"SymmetricMatrixSumTransform"), Transformation(SymmetricMatrixSumTransformParameterNames(n),n), _n(n), _x(), _x_prev(), _any_has_changed(true), _has_changed(), _recalculate(true) {
	}
	
	SymmetricMatrixSumTransform::SymmetricMatrixSumTransform(const SymmetricMatrixSumTransform& x) : DAGcomponent(x), Transformation(x), _n(x._n), _x(x._x), _x_prev(x._x_prev), _any_has_changed(x._any_has_changed), _has_changed(x._has_changed), _recalculate(x._recalculate) {
	}
	
	int SymmetricMatrixSumTransform::length() const {
		if(_recalculate) recalculate();
		return _x.n();
	}
	
	double SymmetricMatrixSumTransform::get_double(const int i, const int j) const {
		if(_recalculate) recalculate();
		return _x.safe(i,j);
	}
	
	vector<double> SymmetricMatrixSumTransform::get_lower_triangle() const {
		if(_recalculate) recalculate();
		vector<double> ret(((length()+1)*length())/2);
		int i,j,k;
		for(i=0,k=0;i<length();i++) {
			for(j=0;j<=i;j++,k++) {
				ret[k] = _x[i][j];
			}
		}
		return ret;
	}
	
	HalfMatrix<double> SymmetricMatrixSumTransform::get_matrix() const {
		if(_recalculate) recalculate();
		return _x;
	}
	
	bool SymmetricMatrixSumTransform::has_changed(const int i, const int j) const {
		if(_recalculate) recalculate();
		return _has_changed.safe(i,j);
	}
	
	bool SymmetricMatrixSumTransform::has_changed() const {
		if(_recalculate) recalculate();
		return _any_has_changed;
	}
	
	bool SymmetricMatrixSumTransform::check_parameter_type(const int i, Variable* parameter) {
		return(dynamic_cast<SymmetricMatrixVariable*>(parameter));
	}
	
	void SymmetricMatrixSumTransform::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
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
	
	void SymmetricMatrixSumTransform::recalculate() const {
		int i,j,k;
		const HalfMatrix<double> &summand = ((const SymmetricMatrixVariable*)get_parameter(0))->get_matrix();
		const int len = summand.n();
		_x_prev = _x;
		_x = summand;
		for(k=1;k<_n;k++) {
			const HalfMatrix<double> &summand = ((const SymmetricMatrixVariable*)get_parameter(k))->get_matrix();
			if(summand.n()!=len) {
				const SymmetricMatrixVariable* v = (const SymmetricMatrixVariable*)get_parameter(k);
				stringstream errMsg;
				errMsg << "SymmetricMatrixSumTransform::recalculate(): SymmetricMatrixVariable object " << v->name() << " has length " << v->length();
				errMsg << " whereas other objects in the sum " << name() << " have length " << len;
				error(errMsg.str().c_str());
			}
			for(i=0;i<len;i++) {
				for(j=0;j<=i;j++) {
					_x[i][j] += summand[i][j];
				}
			}
		}
		// Compare to previous
		if(_x_prev.n()==len) {
			_any_has_changed = false;
			for(i=0;i<len;i++) {
				for(j=0;j<=i;j++) {
					_has_changed[i][j] = _x[i][j]!=_x_prev[i][j];
					if(_has_changed[i][j]) _any_has_changed = true;
				}
			}
		} else {
			_any_has_changed = true;
			_has_changed = HalfMatrix<bool>(len,true);
		}
		_recalculate = false;
	}
	
}	// namespace gcat_lmm

