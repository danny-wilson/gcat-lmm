/*  Copyright 2014 Daniel Wilson.
 *
 *  MultivariateNormal.cpp
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
#include <lmm/Distributions/MultivariateNormal.h>
#include <sstream>
#include <iostream>

using std::stringstream;
using std::cout;
using std::endl;

namespace gcat_lmm {
	
	const string MultivariateNormalDistributionParameterNames[2] = {"mean","cov"};
	
	MultivariateNormalDistribution::MultivariateNormalDistribution(const int length, string name, DAG* dag) : DAGcomponent(name,dag,"MultivariateNormalDistribution"), Distribution(MultivariateNormalDistributionParameterNames,2), _length(length), PI(3.141592653589793238) {
	}
	
	MultivariateNormalDistribution::MultivariateNormalDistribution(const MultivariateNormalDistribution &x) : DAGcomponent(x), Distribution(x), _length(x._length), PI(3.141592653589793238) {
	}
	
	bool MultivariateNormalDistribution::check_random_variable_type(RandomVariable* random_variable) {
		ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(random_variable);
		// Check the length is compatible
		if(!(rv==0) && rv->length()!=length()) {
			stringstream errMsg;
			errMsg << "MultivariateNormalDistribution::check_random_variable_type(): ContinuousVectorVariable object " << rv->name() << " has length " << rv->length();
			errMsg << " whereas MultivariateNormalDistribution object " << name() << " has length " << length();
			error(errMsg.str().c_str());
		}
		return rv;
	}
	
	bool MultivariateNormalDistribution::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0:	{ //	mean
				ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(parameter);
				// Check the length is compatible
				if((!rv==0) && rv->length()!=length()) {
					stringstream errMsg;
					errMsg << "MultivariateNormalDistribution::check_parameter_type(): ContinuousVectorVariable object " << rv->name() << " has length " << rv->length();
					errMsg << " whereas MultivariateNormalDistribution object " << name() << " has length " << length();
					error(errMsg.str().c_str());
				}
				return rv;
			}
			case 1: { //	cov
				SymmetricMatrixVariable* rv = dynamic_cast<SymmetricMatrixVariable*>(parameter);
				// Check the length is compatible
				if((!rv==0) && rv->length()!=length()) {
					stringstream errMsg;
					errMsg << "MultivariateNormalDistribution::check_parameter_type(): SymmetricMatrixVariable object " << rv->name() << " has length " << rv->length();
					errMsg << " whereas MultivariateNormalDistribution object " << name() << " has length " << length();
					error(errMsg.str().c_str());
				}
				return rv;
			}
			default:
				error("MultivariateNormalDistribution::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	int MultivariateNormalDistribution::length() const {
		return _length;
	}	
	
	void MultivariateNormalDistribution::set_mean(ContinuousVectorVariable* mean) {
		set_parameter(0,(Variable*)mean);
	}
	
	void MultivariateNormalDistribution::set_cov(SymmetricMatrixVariable* cov) {
		set_parameter(1,(Variable*)cov);
	}
	
	ContinuousVectorVariable const*  MultivariateNormalDistribution::get_mean() const {
		return (ContinuousVectorVariable const*)get_parameter(0);
	}
	
	SymmetricMatrixVariable const*  MultivariateNormalDistribution::get_cov() const {
		return (SymmetricMatrixVariable const*)get_parameter(1);
	}
	
	mydouble MultivariateNormalDistribution::likelihood(const RandomVariable* rv, const Value* val) {
		if(val==0) error("MultivariateNormalDistribution::likelihood(): variable not found");
		int j,k,l;
		
		// Need functions for computing determinant, inverse and product
		const vector<double> &m = get_mean()->get_doubles();
		//const LowerTriangularMatrix<double> &v = get_cov()->get_lower_triangle();
		const HalfMatrix<double> &v = get_cov()->get_matrix();
		const vector<double> &x = ((ContinuousVectorVariable*)val)->get_doubles();
		const int n = get_mean()->length();
		// Check the lengths are all correct
		if(x.size()!=n) error("MultivariateNormalDistribution::likelihood(): mean did not have same length as random variable");
		if(v.n()!=n) error("MultivariateNormalDistribution::likelihood(): mean did not have same length as covariance matrix");
		// Original matrix
		if(false) {
			cout << "Original matrix:\n";
			for(j=0;j<n;j++) {
				for(k=0;k<n;k++) {
					const double tp = v.safe(j,k);
					cout << tp << "\t";
				}
				cout << endl;
			}
		}
		// Cholesky decomposition
		const bool isposdef = v.Cholesky(_L,false);
		if(false) {
			cout << "Cholesky decomposition:\n";
			for(j=0;j<n;j++) {
				for(k=0;k<n;k++) {
					const double tp = _L.safe(j,k);
					cout << tp << "\t";
				}
				cout << endl;
			}
		}
		
		if(!isposdef) {
			// Temporarily...
			//error("MultivariateNormalDistribution::likelihood(): covariance matrix is not positive definite");
			myutils::warning("MultivariateNormalDistribution::likelihood(): covariance matrix is not positive definite");
			return 0.0;
		}
		
		// Invert lower triangular matrix
		_Linv = _L.invert_lotri();
		if(false) {
			cout << "Inverse Cholesky:\n";
			for(j=0;j<n;j++) {
				for(k=0;k<n;k++) {
					const double tp = _Linv.safe(j,k);
					cout << tp << "\t";
				}
				cout << endl;
			}
		}
		// Compute the log determinant
		const double logdet = log(Cholesky_to_determinant(_L));
		if(false) {
			cout << "Log determinant: " << logdet << endl;
		}
		// Compute exponential term: t(m-x) inv(v) (m-x) = t(m-x) inv(L)' inv(L) (m-x)
		// {t(m-x) inv(L)'}[1,j] = sum_{k=1..n} (m[k]-x[k])*invL[j][k] = sum_{k=1..j} (m[k]-x[k])*invL[j][k]
		// {t(m-x) inv(L)' inv(L)}[1,j] = sum_{l=1..n} {t(m-x) inv(L)'}[1,l] inv(L)[l,j]
		//     = sum_{l=j..n} sum_{k=1..l} (m[k]-x[k])*invL[l][k]*invL[l][j]
		// t(m-x) inv(L)' inv(L) (m-x) = sum_{j=1..n} {t(m-x) inv(L)' inv(L)}[1,j] (m[j]-x[j])
		//     = sum_{j=1..n} sum_{l=j..n} sum_{k=1..l} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j]) ***
		//     = sum_{j=1..n} sum_{k=1..n} sum_{l=1..n} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j])
		//     = sum_{j=1..n} sum_{k=1..n} sum_{l=max(j,k)..n} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j])
		//     = sum_{j=1..n} sum_{k=j..n} sum_{l=k..n} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j])
		double tp = 0.0;
		for(j=0;j<n;j++) {
			for(l=j;l<n;l++) {
				for(k=0;k<=l;k++) {
					tp += (m[k]-x[k])*_Linv[l][k]*_Linv[l][j]*(m[j]-x[j]);
				}
			}
		}
		// pdf: (2*pi)^(-n/2) * det(v)^(-1/2) * exp( -0.5 * t(m-x) %*% inv(v) %*% (m-x))
		mydouble ret;
		ret.setlog(-0.5*((double)n*log(2.0*PI) + logdet + tp));
		if(false) {
			cout << "Log likelihood: " << ret.LOG() << endl;
		}
		return ret;
	}
	
} // namespace gcat_lmm

