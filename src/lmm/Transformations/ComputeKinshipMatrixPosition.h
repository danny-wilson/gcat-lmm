/*  Copyright 2014 Daniel Wilson.
 *
 *  ComputeKinshipMatrix.h
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
#ifndef _COMPUTE_KINSHIP_MATRIX_POSITION_TRANSFORM_H_
#define _COMPUTE_KINSHIP_MATRIX_POSITION_TRANSFORM_H_
#include <DAG/Transformation.h>
#include <Variables/Continuous.h>
#include <Variables/Discrete.h>
#include <Variables/Matrix.h>
#include <lmm/Variables/SymmetricMatrix.h>
#include <halfmatrix.h>

using namespace gcat;
using myutils::HalfMatrix;

namespace gcat_lmm {
	
	class ComputeKinshipMatrixPositionTransform : public SymmetricMatrixVariable, public Transformation {
	private:
		// Initialized?
		mutable bool _init;
		// Recalculate?
		mutable bool _recalculate;
		mutable int _nrecalculate;
		// Internal copy of kinship matrix
		mutable HalfMatrix<double> _x, _x_prev;
		// Keep track of whether elements have changed
		mutable bool _any_has_changed;
		mutable HalfMatrix<bool> _has_changed;
		// Mean frequency per site
		mutable double _f;
	public:
		// Constructor
		ComputeKinshipMatrixPositionTransform(string name="", DAG* dag=0);
		// Copy constructor
		ComputeKinshipMatrixPositionTransform(const ComputeKinshipMatrixPositionTransform& x);
		
		// Implementation of virtual functions inherited from base classes
		int length() const;
		double get_double(const int i, const int j) const;
		vector<double> get_lower_triangle() const;
		HalfMatrix<double> get_matrix() const;
		bool has_changed(const int i, const int j) const;
		bool has_changed() const;
		bool check_parameter_type(const int i, Variable* parameter);
		
		// Convenience functions
		void set_bip(MatrixVariable* vec);
		void set_position(DiscreteVariable* pos);
		MatrixVariable const* get_bip() const;
		DiscreteVariable const* get_position() const;
		
		// Overload method inherited from Transformation
		void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
		
	private:
		void recalculate() const;
	};
	
}	// namespace gcat_lmm

#endif _COMPUTE_KINSHIP_MATRIX_POSITION_TRANSFORM_H_