/*  Copyright 2014 Daniel Wilson.
 *
 *  ComputeKinshipMatrix.cpp
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
#include <lmm/Transformations/ComputeKinshipMatrixPosition.h>
#include <iostream>
#include <ctime>

using std::cout;
using std::endl;
using std::clock;

namespace gcat_lmm {
	
	const string ComputeKinshipMatrixPositionTransformParameterNames[2] = {"bip","position"};
	
	ComputeKinshipMatrixPositionTransform::ComputeKinshipMatrixPositionTransform(string name, DAG* dag) : DAGcomponent(name,dag,"ComputeKinshipMatrixPositionTransform"), Transformation(ComputeKinshipMatrixPositionTransformParameterNames,2), _x(), _x_prev(), _any_has_changed(true), _has_changed(), _recalculate(true), _init(true), _nrecalculate(0) {
	}
	
	ComputeKinshipMatrixPositionTransform::ComputeKinshipMatrixPositionTransform(const ComputeKinshipMatrixPositionTransform& x) : DAGcomponent(x), Transformation(x), _x(x._x), _x_prev(x._x_prev), _any_has_changed(x._any_has_changed), _has_changed(x._has_changed), _recalculate(x._recalculate), _init(x._init), _nrecalculate(x._nrecalculate) {
	}
	
	int ComputeKinshipMatrixPositionTransform::length() const {
		return get_bip()->nrows();
	}
	
	double ComputeKinshipMatrixPositionTransform::get_double(const int i, const int j) const {
		if(_recalculate) recalculate();
		return _x.safe(i,j);
	}
	
	vector<double> ComputeKinshipMatrixPositionTransform::get_lower_triangle() const {
		if(_recalculate) recalculate();
		const int n = length();
		vector<double> ret((n*(n+1))/2);
		int i,j,k;
		for(i=0,k=0;i<n;i++) {
			for(j=0;j<=i;j++,k++) {
				ret[k] = _x[i][j];
			}
		}
		return ret;
	}
	
	HalfMatrix<double> ComputeKinshipMatrixPositionTransform::get_matrix() const {
		if(_recalculate) recalculate();
		return _x;
	}
	
	bool ComputeKinshipMatrixPositionTransform::has_changed(const int i, const int j) const {
		if(_recalculate) recalculate();
		return _any_has_changed;		// Potentially inefficient
	}
	
	bool ComputeKinshipMatrixPositionTransform::has_changed() const {
		if(_recalculate) recalculate();
		return _any_has_changed;
	}
	
	bool ComputeKinshipMatrixPositionTransform::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0:	// Matrix variable
				return(dynamic_cast<MatrixVariable*>(parameter));
			case 1: // Integer
				return(dynamic_cast<DiscreteVariable*>(parameter));
			default:
				error("ComputeKinshipMatrixPositionTransform::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	void ComputeKinshipMatrixPositionTransform::set_bip(MatrixVariable* m) {
		set_parameter(0,(Variable*)m);
	}
	
	void ComputeKinshipMatrixPositionTransform::set_position(DiscreteVariable* pos) {
		set_parameter(1,(Variable*)pos);
	}
	
	MatrixVariable const* ComputeKinshipMatrixPositionTransform::get_bip() const {
		return (MatrixVariable const*)get_parameter(0);
	}
	
	DiscreteVariable const* ComputeKinshipMatrixPositionTransform::get_position() const {
		return (DiscreteVariable const*)get_parameter(1);
	}
	
	void ComputeKinshipMatrixPositionTransform::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
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
	
	void ComputeKinshipMatrixPositionTransform::recalculate() const {
		if(_recalculate) {
			++_nrecalculate;
			std::cout << "ComputeKinshipMatrixPositionTransform::recalculate(): call number " << _nrecalculate << std::endl;
			int i,j;
			const MatrixVariable &bip = *get_bip();
			const DiscreteVariable &position = *get_position();
			const int n = bip.nrows();
			const int L = bip.ncols();
			const int k = position.get_int();
			if(k<0 || k>=bip.ncols()) error("ComputeKinshipMatrixPositionTransform::recalculate(): position outside range");
			_x_prev = _x;
			_x = HalfMatrix<double>(n,0.0);
			time_t start = clock();
			// Calculate mean bip frequency - Bayesian estimate helps (ensure?) the kinship matrix is positive definite
			_f = 1.0;
			for(i=0;i<n;i++) {
				_f += bip.get_double(i,k);
			}
			_f /= (double)(n+2);
			cout << "Calculated bip frequency in " << (clock()-start)/CLOCKS_PER_SEC << " s" << endl;
			start = clock();
			// Do the main computation (slow)
			for(i=0;i<n;i++) {
				const double bip_ik = bip.get_double(i,k);
				for(j=0;j<=i;j++) {
					_x[i][j] += (bip_ik-_f)*(bip.get_double(j,k)-_f);
				}
			}
			cout << "Calculated kinship matrix in " << (clock()-start)/CLOCKS_PER_SEC << " s" << endl;
			/*		For some reason this can break the positive definite property of the kinship matrix 
					Could be a numerical problem. The Bayesian estimate of allele frequency may have got round this.
			 */
			// Divide through by the number of loci
			for(i=0;i<n;i++) {
				for(j=0;j<=i;j++) {
					_x[i][j] /= (double)L;
				}
			}
			if(false) {
				cout << "Kinship matrix:\n";
				const int N = (n>5) ? 5 : n;
				for(i=0;i<N;i++) {
					for(j=0;j<N;j++) {
						cout << _x.safe(i,j) << "\t";
					}
					cout << endl;
				}
			}
			
			// Assume something has changed (inefficient)
			_any_has_changed = true;
			_recalculate = false;
		}
	}
	
}	// namespace gcat_lmm
