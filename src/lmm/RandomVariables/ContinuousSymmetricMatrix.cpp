/*  Copyright 2014 Daniel Wilson.
 *
 *  ContinuousSymmetricMatrix.cpp
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
#include <lmm/RandomVariables/ContinuousSymmetricMatrix.h>

namespace gcat_lmm {

	ContinuousSymmetricMatrixRV::ContinuousSymmetricMatrixRV(const int n, string name, DAG* dag, const HalfMatrix<double>* values) : DAGcomponent(name,dag), RandomVariable(), _n(n), _value(n) {
		if(values!=NULL) {
			if(values->n()!=length()) {
				string errMsg = "ContinuousSymmetricMatrixRV: unexpected number of rows and columns of input matrix for object " + name;
				error(errMsg.c_str());
			}
			int i,j;
			for(i=0;i<_n;i++) {
				for(j=0;j<=i;j++) {
					_value[i][j] = (*values)[i][j];
				}
			}
		}
	}
		
	ContinuousSymmetricMatrixRV::ContinuousSymmetricMatrixRV(const int n, string name, DAG* dag, const vector<double>* x) : DAGcomponent(name,dag), RandomVariable(), _n(n), _value(n) {
		if(x!=NULL) {
			if(x->size()!=((nrows()+1)*ncols())/2) {
				string errMsg = "ContinuousSymmetricMatrixRV: unexpected number of rows/columns of input matrix for object " + name;
				error(errMsg.c_str());
			}
			int i,j,k;
			for(i=0,k=0;i<_n;i++) {
				for(j=0;j<=i;j++,k++) {
					_value[i][j] = (*x)[k];
				}
			}
		}
	}

	ContinuousSymmetricMatrixRV::ContinuousSymmetricMatrixRV(const ContinuousSymmetricMatrixRV& x) : DAGcomponent(x), RandomVariable(x), _n(x._n), _value(x._value) {
	}
	
	ContinuousSymmetricMatrixRV::~ContinuousSymmetricMatrixRV() {};
	
	void ContinuousSymmetricMatrixRV::set(const int i, const int j, const double value) {
		if(j<=i) {
			_value[i][j] = value;
		} else {
			_value[j][i] = value;
		}
	}
	
	void ContinuousSymmetricMatrixRV::set(const HalfMatrix<double>& value) {
		if(value.n()!=length()) {
			string errMsg = "ContinuousSymmetricMatrixRV::set: unexpected number of rows and columns of input matrix for object " + name();
			error(errMsg.c_str());
		}
		int i,j;
		for(i=0;i<nrows();i++) {
			for(j=0;j<=i;j++) {
				_value[i][j] = value[i][j];
			}
		}
	}
	
	int ContinuousSymmetricMatrixRV::length() const {
		return _n;
	}
	
	double ContinuousSymmetricMatrixRV::get_double(const int i, const int j) const {
		if(j<=i) return _value[i][j];
		return _value[j][i];
	}
	
	vector<double> ContinuousSymmetricMatrixRV::get_lower_triangle() const {
		vector<double> ret((_n*(_n+1))/2);
		int i,j,k;
		for(i=0,k=0;i<_n;i++) {
			for(j=0;j<=i;j++,k++) {
				ret[k] = _value[i][j];
			}
		}
		return ret;
	}

	HalfMatrix<double> ContinuousSymmetricMatrixRV::get_matrix() const {
		return _value;
	}

}
