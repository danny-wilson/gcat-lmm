/*  Copyright 2014 Daniel Wilson.
 *
 *  ContinuousSymmetricMatrix.h
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
#ifndef _CONTINUOUS_SYMMETRIC_MATRIX_RANDOM_VARIABLE_H_
#define _CONTINUOUS_SYMMETRIC_MATRIX_RANDOM_VARIABLE_H_
#include <lmm/Variables/SymmetricMatrix.h>
#include <DAG/RandomVariable.h>
#include <halfmatrix.h>

using namespace gcat;
using myutils::HalfMatrix;

namespace gcat_lmm {
	
	/* This class is a Random Variable but currently intended for data, not parameters */
	class ContinuousSymmetricMatrixRV : public SymmetricMatrixVariable, public RandomVariable {
	private:
		// Number of rows and columns
		int _n;
		// Values
		HalfMatrix< double > _value;
	public:
		// Constructor
		ContinuousSymmetricMatrixRV(const int n, string name="", DAG* dag=0, const HalfMatrix<double>* values=NULL);
		// Constructor
		ContinuousSymmetricMatrixRV(const int n, string name="", DAG* dag=0, const vector<double>* x=NULL);
		// Copy constructor
		ContinuousSymmetricMatrixRV(const ContinuousSymmetricMatrixRV& x);
		// Destructor
		virtual ~ContinuousSymmetricMatrixRV();
		
		// Manipulators
		void set(const int i, const int j, const double value);
		void set(const HalfMatrix<double>& value);
//		void set(const vector<int>& posi, const vector<int>& posj, const vector<double>& value);
//		void propose(const int i, const int j, const double value);
//		void propose(const HalfMatrix<double>& value);
//		void propose(const vector<int>& posi, const vector<int>& posj, const vector<double>& value);
//		void accept();
//		void revert();
		
		// Implementation of inherited methods
		int length() const;
		// Get element at particular position
		double get_double(const int i, const int j) const;
		// Get vector of values of the lower triangular matrix
		vector<double> get_lower_triangle() const;
		// Get symmetric matrix object
		HalfMatrix<double> get_matrix() const;
//		bool has_changed(const int i) const;
//		vector<bool> has_changed() const;
		
	};
	
} // namespace gcat_lmm


#endif // _CONTINUOUS_SYMMETRIC_MATRIX_RANDOM_VARIABLE_H_
