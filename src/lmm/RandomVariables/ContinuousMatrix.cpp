/*  Copyright 2014 Daniel Wilson.
 *
 *  ContinuousMatrix.cpp
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
#include <lmm/RandomVariables/ContinuousMatrix.h>

namespace gcat_lmm {

	ContinuousMatrixRV::ContinuousMatrixRV(const int nrow, const int ncol, string name, DAG* dag, const Matrix<double>* values) : DAGcomponent(name,dag), RandomVariable(), _nrows(nrow), _ncols(ncol), _value(nrow,ncol) {
		if(values!=NULL) {
			if(values->nrows()!=nrows()) {
				string errMsg = "ContinuousMatrixRV: unexpected number of rows of input matrix for object " + name;
				error(errMsg.c_str());
			}
			if(values->ncols()!=ncols()) {
				string errMsg = "ContinuousMatrixRV: unexpected number of columns of input matrix for object " + name;
				error(errMsg.c_str());
			}
			int i,j;
			for(i=0;i<_nrows;i++) {
				for(j=0;j<_ncols;j++) {
					_value[i][j] = (*values)[i][j];
				}
			}
		}
	}
		
	ContinuousMatrixRV::ContinuousMatrixRV(const int nrow, const int ncol, string name, DAG* dag, const vector<double>* x) : DAGcomponent(name,dag), RandomVariable(), _nrows(nrow), _ncols(ncol), _value(nrow,ncol) {
		if(x!=NULL) {
			if(x->size()!=(nrows()*ncols())) {
				string errMsg = "ContinuousMatrixRV: unexpected number of rows/columns of input matrix for object " + name;
				error(errMsg.c_str());
			}
			int i,j,k;
			for(i=0,k=0;i<_nrows;i++) {
				for(j=0;j<_ncols;j++,k++) {
					_value[i][j] = (*x)[k];
				}
			}
		}
	}

	ContinuousMatrixRV::ContinuousMatrixRV(const ContinuousMatrixRV& x) : DAGcomponent(x), RandomVariable(x), _nrows(x._nrows), _ncols(x._ncols), _value(x._value) {
	}
	
	ContinuousMatrixRV::~ContinuousMatrixRV() {};
	
	void ContinuousMatrixRV::set(const int i, const int j, const double value) {
		_value[i][j] = value;
	}
	
	void ContinuousMatrixRV::set(const Matrix<double>& value) {
		if(value.nrows()!=nrows()) {
			string errMsg = "ContinuousMatrixRV::set: unexpected number of rows of input matrix for object " + name();
			error(errMsg.c_str());
		}
		if(value.ncols()!=ncols()) {
			string errMsg = "ContinuousMatrixRV::set: unexpected number of columns of input matrix for object " + name();
			error(errMsg.c_str());
		}
		int i,j;
		for(i=0;i<nrows();i++) {
			for(j=0;j<ncols();j++) {
				_value[i][j] = value[i][j];
			}
		}
	}
	
	int ContinuousMatrixRV::nrows() const {
		return _nrows;
	}
	
	int ContinuousMatrixRV::ncols() const {
		return _ncols;
	}
	
	double ContinuousMatrixRV::get_double(const int i, const int j) const {
		return _value[i][j];
	}

}
