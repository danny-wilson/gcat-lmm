/*  Copyright 2014 Daniel Wilson.
 *
 *  RescaleVector.h
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
#ifndef _RESCALE_VECTOR_TRANSFORM_H_
#define _RESCALE_VECTOR_TRANSFORM_H_
#include <DAG/Transformation.h>
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>

using namespace gcat;

namespace gcat_lmm {
	
	class RescaleVectorTransform : public ContinuousVectorVariable, public Transformation {
	private:
		// Initialized?
		mutable bool _init;
		// Recalculate?
		mutable bool _recalculate;
		// Internal copy of concatenated items
		mutable vector<double> _x, _x_prev;
		// Keep track of whether elements have changed
		mutable vector<bool> _has_changed;
	public:
		// Constructor
		RescaleVectorTransform(string name="", DAG* dag=0);
		// Copy constructor
		RescaleVectorTransform(const RescaleVectorTransform& x);
		
		// Implementation of virtual functions inherited from base classes
		int length() const;
		double get_double(const int i) const;
		vector<double> get_doubles() const;
		bool has_changed(const int i) const;
		vector<bool> has_changed() const;
		bool check_parameter_type(const int i, Variable* parameter);
		
		// Convenience functions
		void set_scalar(ContinuousVariable* sca);
		void set_vector(ContinuousVectorVariable* vec);
		ContinuousVariable const* get_scalar() const;
		ContinuousVectorVariable const* get_vector() const;

		// Overload method inherited from Transformation
		void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
		
	private:
		void recalculate() const;
	};
	
}	// namespace gcat_lmm

#endif _RESCALE_VECTOR_TRANSFORM_H_