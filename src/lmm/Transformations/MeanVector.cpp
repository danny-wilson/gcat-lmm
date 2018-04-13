/*  Copyright 2014 Daniel Wilson.
 *
 *  MeanVector.cpp
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
#include <lmm/Transformations/MeanVector.h>

namespace gcat_lmm {
	
	const string MeanVectorParameterNames[1] = {"x"};
	
	MeanVector::MeanVector(string name, DAG* dag) : DAGcomponent(name,dag,"MeanVector"), Transformation(MeanVectorParameterNames,1) {
	}
	
	MeanVector::MeanVector(const MeanVector& x) : DAGcomponent(x), Transformation(x) {
	}
	
	double MeanVector::get_double() const {
		const int n = get_vector()->length();
		double sum = 0.0;
		int i;
		for(i=0;i<n;i++) {
			sum += get_vector()->get_double(i);
		}
		return sum/(double)n;
	}
	
	bool MeanVector::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0:
				return(dynamic_cast<ContinuousVectorVariable*>(parameter));
			default:
				error("MeanVector::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	void MeanVector::set_vector(ContinuousVectorVariable* v) {
		set_parameter(0,(Variable*)v);
	}
	
	ContinuousVectorVariable const* MeanVector::get_vector() const {
		return (ContinuousVectorVariable const*)get_parameter(0);
	}
	
} // namespace gcat_lmm

