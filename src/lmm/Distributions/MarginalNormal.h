/*  Copyright 2014 Daniel Wilson.
 *
 *  MarginalNormal.h
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
#ifndef _MARGINAL_NORMAL_DISTRIBUTION_H_
#define _MARGINAL_NORMAL_DISTRIBUTION_H_
#include <DAG/Distribution.h>
#include <RandomVariables/Continuous.h>
#include <RandomVariables/ContinuousVector.h>
#include <lmm/Variables/SymmetricMatrix.h>
#include <Properties/Length.h>
#include <halfmatrix.h>

using namespace gcat;
using myutils::HalfMatrix;

namespace gcat_lmm {
	
	class MarginalNormalDistribution : public Distribution, public LengthProperty {
	private:
		const double logPI;
		const int _length;
		HalfMatrix<double> _L, _Linv;
	public:
		// Constructor
		MarginalNormalDistribution(const int length, string name="", DAG* dag=0);
		// Copy constructor
		MarginalNormalDistribution(const MarginalNormalDistribution &x);
		// Implementation of virtual function inherited from base class Distribution
		bool check_random_variable_type(RandomVariable* random_variable);
		// Implementation of virtual function inherited from class DependentVariable
		bool check_parameter_type(const int i, Variable* parameter);
		// Implementation of inherited function
		int length() const;
		// Convenience functions
		void set_cov(SymmetricMatrixVariable* cov);
		SymmetricMatrixVariable const* get_cov() const;
		
		// Compute log-likelihood
		mydouble likelihood(const RandomVariable* rv, const Value* val);
	};
	
} // namespace gcat_lmm


#endif //_MARGINAL_NORMAL_DISTRIBUTION_H_


