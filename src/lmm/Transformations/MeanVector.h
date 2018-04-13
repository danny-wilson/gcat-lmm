/*  Copyright 2014 Daniel Wilson.
 *
 *  MeanVector.h
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
#ifndef _MEAN_VECTOR_H_
#define _MEAN_VECTOR_H_
#include <Variables/Continuous.h>
#include <DAG/Transformation.h>
#include <Variables/ContinuousVector.h>

using namespace gcat;

namespace gcat_lmm {
	
	class MeanVector : public ContinuousVariable, public Transformation {
	public:
		// Constructor
		MeanVector(string name="", DAG* dag=0);
		// Copy constructor
		MeanVector(const MeanVector& x);
		
		// Implementation of virtual functions inherited from base classes
		double get_double() const;
		bool check_parameter_type(const int i, Variable* parameter);
		
		// Convenience functions
		void set_vector(ContinuousVectorVariable* v);
		ContinuousVectorVariable const* get_vector() const;	
	};
	
} // namespace gcat_lmm

#endif // _MEAN_VECTOR_H_
