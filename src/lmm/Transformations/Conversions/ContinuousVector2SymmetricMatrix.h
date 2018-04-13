/*  Copyright 2014 Daniel Wilson.
 *
 *  ContinuousVector2SymmetricMatrix.h
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
#ifndef _CONTINUOUS_VECTOR_VARIABLE_2_SYMMETRIC_MATRIX_VARIABLE_H_
#define _CONTINUOUS_VECTOR_VARIABLE_2_SYMMETRIC_MATRIX_VARIABLE_H_
#include <Variables/ContinuousVector.h>
#include <lmm/Variables/SymmetricMatrix.h>
#include <DAG/Transformation.h>

namespace gcat_lmm {
	
	const string ContinuousVectorVariable2SymmetricMatrixVariableParameterNames[1] = {"x"};
	
	class ContinuousVectorVariable2SymmetricMatrixVariable : public SymmetricMatrixVariable, public Transformation {
	public:
		// Constructor
		ContinuousVectorVariable2SymmetricMatrixVariable(string name="", DAG* dag=0) : DAGcomponent(name,dag,"ContinuousVectorVariable2SymmetricMatrixVariable"), Transformation(ContinuousVectorVariable2SymmetricMatrixVariableParameterNames,1) {};
		// Copy constructor
		ContinuousVectorVariable2SymmetricMatrixVariable(const ContinuousVectorVariable2SymmetricMatrixVariable& x) : DAGcomponent(x), Transformation(x) {};
		
		// Implementation of virtual functions inherited from MatrixVariable
		// Get value
		double get_double(int i, int j) const {			
			if(j>i) {
				// Ensure j<=i
				SWAP(i,j);
			}
			// First element of the row follows an arithmetic sum starting at 0 with difference 1
			// 0 (+1), 1 (+2), 3 (+3), 6 (+4), 10 (+5), 15 (+6), 21
			// Given by equation (i+1)*((i+1)-1)/2, where i is the (0-based) row
			// So k = (i+1)*i/2 + j, for j<=i
			// Therefore
			const int k = ((i+1)*i)/2 + j;
			return get_x()->get_double(k);
		}
		
		// Implementation of virtual functions inherited from SymmetricMatrixVariable
		// Get length of the variable
		int length() const {
			// Lower triangle has length l = n*(n+1)/2
			// So n = (sqrt(8*l+1)-1)/2
			const double l = (double)get_x()->length();
			const int n = (int)round(sqrt(8.0*l+1.0)-1)/2;
			return n;
		}		
		// Get vector of values of the lower triangular matrix
		vector<double> get_lower_triangle() const {
			return get_x()->get_doubles();
		}
		// Get symmetric matrix object
		HalfMatrix<double> get_matrix() const {
			const int n = length();
			HalfMatrix<double> ret(n,myutils::SYMMETRIC);
			int i,j,k;
			for(i=0,k=0;i<n;i++) {
				for(j=0;j<=i;j++,k++) {
					ret[i][j] = get_x()->get_double(k);
				}
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
					return(dynamic_cast<ContinuousVectorVariable*>(parameter));
				default:
					error("ContinuousVectorVariable2SymmetricMatrixVariable::check_parameter_type(): parameter not found");
			}
			return false;
		}
		
		// Convenience functions
		void set_x(ContinuousVectorVariable* x) {
			set_parameter(0,(Variable*)x);
		}
		ContinuousVectorVariable const* get_x() const {
			return (ContinuousVectorVariable const*)get_parameter(0);
		}
	};

} // namespace gcat

#endif // _CONTINUOUS_VECTOR_VARIABLE_2_SYMMETRIC_MATRIX_VARIABLE_H_
