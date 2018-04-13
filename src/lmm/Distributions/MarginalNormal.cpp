/*  Copyright 2014 Daniel Wilson.
 *
 *  MarginalNormal.cpp
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
#include <lmm/Distributions/MarginalNormal.h>
#include <sstream>
#include <iostream>
#include <cmath>

using std::stringstream;
using std::cout;
using std::endl;

namespace gcat_lmm {
	
	const string MarginalNormalDistributionParameterNames[1] = {"cov"};
	
	MarginalNormalDistribution::MarginalNormalDistribution(const int length, string name, DAG* dag) : DAGcomponent(name,dag,"MarginalNormalDistribution"), Distribution(MarginalNormalDistributionParameterNames,1), _length(length), logPI(1.1447298858494) {
	}
	
	MarginalNormalDistribution::MarginalNormalDistribution(const MarginalNormalDistribution &x) : DAGcomponent(x), Distribution(x), _length(x._length), logPI(1.1447298858494) {
	}
	
	bool MarginalNormalDistribution::check_random_variable_type(RandomVariable* random_variable) {
		ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(random_variable);
		// Check the length is compatible
		if(!(rv==0) && rv->length()!=length()) {
			stringstream errMsg;
			errMsg << "MarginalNormalDistribution::check_random_variable_type(): ContinuousVectorVariable object " << rv->name() << " has length " << rv->length();
			errMsg << " whereas MarginalNormalDistribution object " << name() << " has length " << length();
			error(errMsg.str().c_str());
		}
		return rv;
	}
	
	bool MarginalNormalDistribution::check_parameter_type(const int i, Variable* parameter) {
		switch(i) {
			case 0: { //	cov
				SymmetricMatrixVariable* rv = dynamic_cast<SymmetricMatrixVariable*>(parameter);
				// Check the length is compatible
				if((!rv==0) && rv->length()!=length()) {
					stringstream errMsg;
					errMsg << "MarginalNormalDistribution::check_parameter_type(): SymmetricMatrixVariable object " << rv->name() << " has length " << rv->length();
					errMsg << " whereas MarginalNormalDistribution object " << name() << " has length " << length();
					error(errMsg.str().c_str());
				}
				return rv;
			}
			default:
				error("MarginalNormalDistribution::check_parameter_type(): parameter not found");
		}
		return false;
	}
	
	int MarginalNormalDistribution::length() const {
		return _length;
	}	
	
	void MarginalNormalDistribution::set_cov(SymmetricMatrixVariable* cov) {
		set_parameter(0,(Variable*)cov);
	}
	
	SymmetricMatrixVariable const*  MarginalNormalDistribution::get_cov() const {
		return (SymmetricMatrixVariable const*)get_parameter(0);
	}
	
	mydouble MarginalNormalDistribution::likelihood(const RandomVariable* rv, const Value* val) {
		if(val==0) error("MarginalNormalDistribution::likelihood(): variable not found");
		int j,k,l;
		
		// Get the parameters and random variable
		const HalfMatrix<double> &V = get_cov()->get_matrix();
		const int n = V.n();
		vector<double> y = ((ContinuousVectorVariable*)val)->get_doubles();
		// Zero-centre y
		double ybar = 0.0;
		for(j=0;j<n;j++) ybar += y[j];
		ybar /= (double)n;
		for(j=0;j<n;j++) y[j] -= ybar;
		// Check the lengths are all correct
		if(y.size()!=n) error("MarginalNormalDistribution::likelihood(): covariance matrix did not have same dimension as random variable");
		// Original matrix
		if(false) {
			cout << "Original matrix:\n";
			const int N = (n>5) ? 5 : n;
			for(j=0;j<N;j++) {
				for(k=0;k<N;k++) {
					const double tp = V.safe(j,k);
					cout << tp << "\t";
				}
				cout << endl;
			}
		}
		// Cholesky decomposition
		const bool isposdef = V.Cholesky(_L,false);
		if(false) {
			cout << "Cholesky decomposition:\n";
			const int N = (n>5) ? 5 : n;
			for(j=0;j<N;j++) {
				for(k=0;k<N;k++) {
					const double tp = _L.safe(j,k);
					cout << tp << "\t";
				}
				cout << endl;
			}
		}
		
		if(!isposdef) {
			// Temporarily...
			//error("MarginalNormalDistribution::likelihood(): covariance matrix is not positive definite");
			myutils::warning("MarginalNormalDistribution::likelihood(): covariance matrix is not positive definite");
			return 0.0;
		}
		
		// Invert lower triangular matrix
		_Linv = _L.invert_lotri();
		if(false) {
			cout << "Inverse Cholesky:\n";
			const int N = (n>5) ? 5 : n;
			for(j=0;j<N;j++) {
				for(k=0;k<N;k++) {
					const double tp = _Linv.safe(j,k);
					cout << tp << "\t";
				}
				cout << endl;
			}
		}
		// Compute the log determinant of V
		const double logdet = log(Cholesky_to_determinant(_L));
		if(false) {
			cout << "Log determinant: " << logdet << endl;
		}
		// Compute "exponential" term: t(m-x) inv(v) (m-x) = t(m-x) inv(L)' inv(L) (m-x)
		// (as per MultivariateNormal.cpp)
		// {t(m-x) inv(L)'}[1,j] = sum_{k=1..n} (m[k]-x[k])*invL[j][k] = sum_{k=1..j} (m[k]-x[k])*invL[j][k]
		// {t(m-x) inv(L)' inv(L)}[1,j] = sum_{l=1..n} {t(m-x) inv(L)'}[1,l] inv(L)[l,j]
		//     = sum_{l=j..n} sum_{k=1..l} (m[k]-x[k])*invL[l][k]*invL[l][j]
		// t(m-x) inv(L)' inv(L) (m-x) = sum_{j=1..n} {t(m-x) inv(L)' inv(L)}[1,j] (m[j]-x[j])
		//     = sum_{j=1..n} sum_{l=j..n} sum_{k=1..l} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j]) ***
		//     = sum_{j=1..n} sum_{k=1..n} sum_{l=1..n} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j])
		//     = sum_{j=1..n} sum_{k=1..n} sum_{l=max(j,k)..n} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j])
		//     = sum_{j=1..n} sum_{k=j..n} sum_{l=k..n} (m[k]-x[k])*invL[l][k]*invL[l][j]*(m[j]-x[j])
		double Q = 0.0;
		for(j=0;j<n;j++) {
			for(l=j;l<n;l++) {
				for(k=0;k<=l;k++) {
					Q += y[k]*_Linv[l][k]*_Linv[l][j]*y[j];
				}
			}
		}
		// Compute the marginal likelihood (up to a constant because we're assuming an improper prior)
		// pdf: lgamma(p/2) -1/2*determinant(V,log=TRUE)$modulus -p/2*log(pi) -p/2*log(Q)
		mydouble ret;
		const double p = (double)n;
		ret.setlog(lgamma(0.5*p) -0.5*logdet -0.5*p*logPI -0.5*p*log(Q));
		if(false) {
			cout << "Log likelihood: " << ret.LOG() << endl;
			cout << lgamma(0.5*p) << "\t" << -0.5*logdet << "\t" << -0.5*p*logPI << "\t" << -0.5*p*log(Q) << endl;
		}
		return ret;
	}
	
} // namespace gcat_lmm

