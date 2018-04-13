/*  Copyright 2014 Daniel Wilson.
 *
 *  ContinuousVariable2DiagonalMatrix.h
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
#ifndef _CONTINUOUS_VARIABLE_2_DIAGONAL_MATRIX_VARIABLE_H_
#define _CONTINUOUS_VARIABLE_2_DIAGONAL_MATRIX_VARIABLE_H_
#include <Variables/Continuous.h>
#include <lmm/Variables/SymmetricMatrix.h>
#include <DAG/Transformation.h>

namespace gcat_lmm {
	
	const string ContinuousVariable2DiagonalMatrixVariableParameterNames[1] = {"x"};
	
	class ContinuousVariable2DiagonalMatrixVariable : public SymmetricMatrixVariable, public Transformation {
	private:
		const int _length;
	public:
		// Constructor
		ContinuousVariable2DiagonalMatrixVariable(const int length, string name="", DAG* dag=0) : DAGcomponent(name,dag,"ContinuousVariable2DiagonalMatrixVariable"), Transformation(ContinuousVariable2DiagonalMatrixVariableParameterNames,1), _length(length) {};
		// Copy constructor
		ContinuousVariable2DiagonalMatrixVariable(const ContinuousVariable2DiagonalMatrixVariable& x) : DAGcomponent(x), Transformation(x), _length(x._length) {};
		
		// Implementation of virtual functions inherited from MatrixVariable
		// Get value
		double get_double(int i, int j) const {
			if(i!=j) return 0.0;
			return get_x()->get_double();
		}
		
		// Implementation of virtual functions inherited from SymmetricMatrixVariable
		// Get length of the variable
		int length() const {
			return _length;
		}		
		// Get vector of values of the lower triangular matrix
		vector<double> get_lower_triangle() const {
			vector<double> ret(((_length+1)*_length)/2,0.0);
			const double diag = get_x()->get_double();
			int i,j,k;
			for(i=0,k=0;i<_length;i++) {
				for(j=0;j<=i;j++,k++) {
					if(i==j) {
						ret[k] = diag;
					}
				}
			}
			return ret;
		}
		// Get symmetric matrix object
		HalfMatrix<double> get_matrix() const {
			const int n = length();
			HalfMatrix<double> ret(n,0.0,myutils::SYMMETRIC);
			const double diag = get_x()->get_double();
			int i;
			for(i=0;i<n;i++) {
				ret[i][i] = diag;
			}
			return ret;
		}
		// Has the value changed at all?
		bool has_changed() const {
			return true;	// Inefficient!
		}
		// Has the value changed at position i, j?
		bool has_changed(const int i, const int j) const {
			return true;	// Inefficient!
		}
			
		bool check_parameter_type(const int i, Variable* parameter)  {
			switch(i) {
				case 0:	// x
					return(dynamic_cast<ContinuousVariable*>(parameter));
				default:
					error("ContinuousVariable2DiagonalMatrixVariable::check_parameter_type(): parameter not found");
			}
			return false;
		}
		
		// Convenience functions
		void set_x(ContinuousVariable* x) {
			set_parameter(0,(Variable*)x);
		}
		ContinuousVariable const* get_x() const {
			return (ContinuousVariable const*)get_parameter(0);
		}
	};

} // namespace gcat

#endif // _CONTINUOUS_VARIABLE_2_DIAGONAL_MATRIX_VARIABLE_H_
