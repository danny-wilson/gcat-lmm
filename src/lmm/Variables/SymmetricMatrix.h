/*  Copyright 2014 Daniel Wilson.
 *
 *  SymmetricMatrix.h
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
// NB:- assumed Continuous
#ifndef _SYMMETRIC_MATRIX_VARIABLE_H_
#define _SYMMETRIC_MATRIX_VARIABLE_H_
#include <Variables/Matrix.h>
#include <myerror.h>
#include <string>
#include <ostream>
#include <halfmatrix.h>
#include <Properties/Length.h>
#include <vector>

using std::string;
using std::ostream;
using std::vector;

using namespace gcat;
using myutils::HalfMatrix;

namespace gcat_lmm {
	
	// Abstract base class derived from MatrixVariable
	// Guarantees methods called get_double(const int i, const int j), etc and implements print() methods
	class SymmetricMatrixVariable : public MatrixVariable, public LengthProperty {
	public:
		// Constructor
		SymmetricMatrixVariable() {};
		// Copy constructor
		SymmetricMatrixVariable(const SymmetricMatrixVariable &x) {};
		// Destructor
		virtual ~SymmetricMatrixVariable() {};
		// Check if the matrix is positive definite as promised (with default method)
		//virtual bool confirm_posdef() const {
			// Failure to confirm does not mean it is not positive definite
		//	return false;
		//}
		// Implementation of number of rows
		virtual int nrows() const {
			return length();
		}
		// Implementation of number of columns
		virtual int ncols() const {
			return length();
		}
		// Get vector of values of the lower triangular matrix
		virtual vector<double> get_lower_triangle() const = 0;
		// Get symmetric matrix object
		virtual HalfMatrix<double> get_matrix() const = 0;
		// Has the value changed at all?
		virtual bool has_changed() const = 0;
		// Has the value changed at position i, j?
		virtual bool has_changed(const int i, const int j) const = 0;
		// Print header (implementation of inherited method)
		virtual void print_header(ostream& out, string sep) {
			int i,j;
			for(i=0;i<nrows();i++) {
				for(j=0;j<=i;j++) {
					if(i>0 || j>0) out << sep;
					out << name() << i << ":" << j;
				}
			}
		}
		// Print value (implementation of inherited method)
		virtual void print(ostream& out, string sep) {
			int i,j;
			for(i=0;i<nrows();i++) {
				for(j=0;j<=i;j++) {
					if(i>0 || j>0) out << sep;
					try{
						out << get_double(i,j);
					}
					catch (BadValueException &e) {
						out << "NA";
					}
				}
			}
		}
	};
	
} // namespace gcat_lmm

#endif // _SYMMETRIC_MATRIX_VARIABLE_H_

