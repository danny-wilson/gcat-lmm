/*  Copyright 2014 Daniel Wilson.
 *
 *  PositiveDefiniteMatrix.h
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
#include <lmm/lmmXML.h>
#include <Distributions/DistributionsXML.h>
#include <RandomVariables/RandomVariablesXML.h>
#include <Transformations/TransformationsXML.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>
#include <lmm/gcat_lmm1.0.xsd.h>
#include <Properties/Length.h>
#include <lmm/Transformations/Conversions/ContinuousVector2SymmetricMatrix.h>
#include <lmm/Distributions/MultivariateNormal.h>
#include <string>
#include <lmm/RandomVariables/ContinuousMatrix.h>
#include <Variables/Continuous.h>
#include <lmm/Transformations/RescaleVector.h>
#include <lmm/Transformations/RescaleSymmetricMatrix.h>
#include <lmm/Transformations/ComputeKinshipMatrix.h>
#include <lmm/Transformations/ComputeKinshipMatrixPosition.h>
#include <lmm/Transformations/Conversions/ContinuousVariable2DiagonalMatrix.h>
#include <RandomVariables/Continuous.h>
#include <RandomVariables/Discrete.h>
#include <lmm/Transformations/SymmetricMatrixSumTransform.h>
#include <lmm/Transformations/Centre.h>
#include <lmm/Transformations/MeanVector.h>
#include <lmm/Distributions/MarginalNormal.h>
#include <ctime>
#include <lmm/Distributions/ContinuousVectorMixture.h>
#include <lmm/Transformations/ContinuousVectorMixtureComponentLogLikelihood.h>

using namespace gcat;
using std::clock;

namespace gcat_lmm {
	
	// DISTRIBUTIONS
	continuous_vector_mixture_distribution_XMLParser::continuous_vector_mixture_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<continuous_vector_mixture_distribution_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 3;
		const char* attrNames[nattr] = {"id","length","p"};
		sattr = attributesToStrings(nattr,attrNames,attrs);
		// Don't instantiate the variable until the values have been read in
	}
	
	void continuous_vector_mixture_distribution_XMLParser::implement_characters(const XMLCh* const chars, const XMLSize_t length) {
		string message = "";
		int i;
		for(i=0;i<length;i++) message += chars[i];
		vector<string> distributions;
		if(!string_to_vector(distributions,message)) error("continuous_vector_mixture_distribution_XMLParser: error interpretting input");
		if(distributions.size()==0) error("continuous_vector_mixture_distribution_XMLParser: no distributions entered");
		const int n = distributions.size();
		// Convert length to int
		int int_length;
		if(!from_string<int>(int_length,sattr[1])) {
			Parameter* rv = getDAG()->get_parameter(sattr[1]);
			if(rv==0) error("continuous_vector_mixture_distribution_XMLParser: could not convert length to int nor find named variable");
			LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
			if(lp==0) error("continuous_vector_mixture_distribution_XMLParser: named variable does not have length property");
			int_length = lp->length();
		}
		if(int_length<=0) error("continuous_vector_mixture_distribution_XMLParser: length must be a positive integer");
		new ContinuousVectorMixture(int_length,n,sattr[0],getDAG());
		for(i=0;i<n;i++) {
			stringstream s;
			s << "distribution" << i;
			getDAG()->assign_distribution_to_compound_distribution(sattr[0],s.str(),distributions[i]);
		}
		// Create the parameter if the default is specified
		if(sattr[2]=="") {
			stringstream name;
			name << "_" << sattr[0];
			const vector<double> p_default(n,1.0/(double)n);
			new ContinuousVectorRV(n,name.str(),getDAG(),p_default);
			getDAG()->set_constant(name.str());
			sattr[2] = name.str();
		}
		getDAG()->assign_parameter_to_distribution(sattr[0],"p",sattr[2]);
	}
	
	marginal_normal_distribution_XMLParser::marginal_normal_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<marginal_normal_distribution_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 4;
		const char* attrNames[nattr] = {"id","distribution","length","cov"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Convert length to int
		int int_length;
		if(!from_string<int>(int_length,sattr[2])) {
			Parameter* rv = getDAG()->get_parameter(sattr[2]);
			if(rv==0) error("marginal_normal_distribution_XMLParser: could not convert length to int nor find named variable");
			LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
			if(lp==0) error("marginal_normal_distribution_XMLParser: named variable does not have length property");
			int_length = lp->length();
		}
		if(int_length<=0) error("marginal_normal_distribution_XMLParser: length must be a positive integer");
		// Instantiate the variable
		new MarginalNormalDistribution(int_length,sattr[0],getDAG());
		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[3],sattr[3]);
	}
	
	multivariate_normal_distribution_XMLParser::multivariate_normal_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<multivariate_normal_distribution_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 5;
		const char* attrNames[nattr] = {"id","distribution","length","mean","cov"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Convert length to int
		int int_length;
		if(!from_string<int>(int_length,sattr[2])) {
			Parameter* rv = getDAG()->get_parameter(sattr[2]);
			if(rv==0) error("multivariate_normal_distribution_XMLParser: could not convert length to int nor find named variable");
			LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
			if(lp==0) error("multivariate_normal_distribution_XMLParser: named variable does not have length property");
			int_length = lp->length();
		}
		if(int_length<=0) error("multivariate_normal_distribution_XMLParser: length must be a positive integer");
		new MultivariateNormalDistribution(int_length,sattr[0],getDAG());
		int i;
		for(i=3;i<nattr;i++) {
			getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
		}
	}
	
	// RANDOM VARIABLES
	continuous_vector_file_XMLParser::continuous_vector_file_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser, const bool constant) : DAGXMLParserTemplate<continuous_vector_file_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 4;
		const char* attrNames[nattr] = {"id","distribution","file","sep"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Read values from file
		vector<double> val(0);
		ifstream in(sattr[2].c_str());
		if(!in.good()) {
			string errMsg = "continuous_vector_file_XMLParser: file " + sattr[2] + " not found";
			error(errMsg.c_str());
		}
		if(sattr[3]!="") error("continuous_vector_file_XMLParser: non-default delimiters not yet implemented");
		double dblin;
		in >> dblin;
		while(!in.eof()) {
			if(in.bad()) {
				string errMsg = "continuous_vector_file_XMLParser: could not read decimal in " + sattr[2];
				error(errMsg.c_str());
			}
			if(!in.good()) {
				string errMsg = "continuous_vector_file_XMLParser: unexpected problem reading " + sattr[2];
				error(errMsg.c_str());
			}
			val.push_back(dblin);
			in >> dblin;
		}
		// Instantiate the variable
		new ContinuousVectorRV(val.size(),sattr[0],getDAG(),val);
		if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],attrNames[1],sattr[1]);
		if(constant) getDAG()->set_constant(sattr[0]);
	}
	
	continuous_matrix_file_XMLParser::continuous_matrix_file_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser, const bool constant) : DAGXMLParserTemplate<continuous_matrix_file_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 4;
		const char* attrNames[nattr] = {"id","distribution","file","sep"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Read values from file
		time_t start = clock();
		vector<double> val(0);
		ifstream in(sattr[2].c_str());
		if(!in.good()) {
			string errMsg = "continuous_matrix_file_XMLParser: file " + sattr[2] + " not found";
			error(errMsg.c_str());
		}
		if(sattr[3]!="") error("continuous_matrix_file_XMLParser: non-default delimiters not yet implemented");
		// Read the file one line at a time: impose constraint that the number of elements is the same for all lines
		double dblin;
		string line;
		getline(in,line);
		bool firstline = true;
		int nrows=0, ncols;
		while(!in.eof()) {
			if(!in.good()) {
				string errMsg = "continuous_matrix_file_XMLParser: unexpected problem reading " + sattr[2];
				error(errMsg.c_str());
			}
			// Convert to vector of doubles
			istringstream linein(line);
			int n = 0;
			while(linein >> dblin) {
				if(linein.bad()) {
					string errMsg = "continuous_matrix_file_XMLParser: could not read decimal in " + sattr[2];
					error(errMsg.c_str());
				}
				val.push_back(dblin);
				++n;
			}
			if(firstline) {
				ncols = n;
				firstline = false;
			} else {
				if(ncols!=n) {
					stringstream errMsg;
					errMsg << "continuous_matrix_file_XMLParser: nrows " << nrows+1 << " contained " << n << " columns, not " << ncols << " as expected";
					error(errMsg.str().c_str());
				}
			}
			++nrows;
			getline(in,line);
		}
		cout << "Read " << sattr[2] << " in " << (clock()-start)/CLOCKS_PER_SEC << " s" << endl;
		// Instantiate the variable
		if(val.size()!=(nrows*ncols)) error("continuous_matrix_file_XMLParser: # values read does not match # rows and columns");
		new ContinuousMatrixRV(nrows,ncols,sattr[0],getDAG(),&val);
		if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],attrNames[1],sattr[1]);
		if(constant) getDAG()->set_constant(sattr[0]);
	}
	
	// TRANSFORMATIONS: CONVERSIONS

	to_diagonal_matrix_XMLParser::to_diagonal_matrix_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<to_diagonal_matrix_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 3;
		const char* attrNames[nattr] = {"id","x","length"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// If length is blank, by default it is interpreted to mean the length of x
		if(sattr[2]=="") sattr[2] = sattr[1];
		// Allow the length to be specified numerically, otherwise obtain it as a property from another object
		int int_length;
		if(!from_string<int>(int_length,sattr[2])) {
			Parameter* rv = getDAG()->get_parameter(sattr[2]);
			if(rv==0) error("to_diagonal_matrix_XMLParser: could not convert length to int nor find named variable");
			LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
			if(lp==0) error("to_diagonal_matrix_XMLParser: named variable does not have length property");
			int_length = lp->length();
		}
		if(int_length<=0) error("to_diagonal_matrix_XMLParser: length must be a positive integer");
		// Allow the variable under transformation to be specified numerically in the "x" element
		double double_val;
		if(from_string(double_val,sattr[1])) {
			// Create a dummy variable with an internally-generated name. Since it won't be possible to refer to this variable by name, force it to be constant
			// It will not be possible to put a distribution on a constant created this way
			stringstream name;
			name << "_" << sattr[0];
			new ContinuousRV(name.str(),getDAG(),double_val);
			getDAG()->set_constant(name.str());
			
			new ContinuousVariable2DiagonalMatrixVariable(int_length,sattr[0],getDAG());
			getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],name.str());
			
			return;
		}
		// Alternatively allow x to specify a named ContinuousVariable
		// Obtain random variable
		Variable* v = getDAG()->get_variable(sattr[1]);
		// Try to convert from ContinuousVectorVariable
		ContinuousVariable* cv = dynamic_cast<ContinuousVariable*>(v);
		if(cv!=0) {
			new ContinuousVariable2DiagonalMatrixVariable(int_length,sattr[0],getDAG());
			getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
			return;
		}
		// Unimplemented conversions: ContinuousVectorVariable		
		// Try other conversions here, in decreasing order of preference
		string errMsg = "to_diagonal_matrix_XMLParser: could not convert variable " + sattr[1] + " of type " + v->type() + " to diagonal matrix";
		error(errMsg.c_str());
	}

	to_symmetric_matrix_XMLParser::to_symmetric_matrix_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<to_symmetric_matrix_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 2;
		const char* attrNames[nattr] = {"id","x"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Obtain random variable
		Variable* v = getDAG()->get_variable(sattr[1]);
		// Test if already of desired type
		if(dynamic_cast<SymmetricMatrixVariable*>(v)!=0) {
			string errMsg = "to_symmetric_matrix_XMLParser: variable " + sattr[1] + " is already of desired type";
			error(errMsg.c_str());
		}
		// Try to convert from ContinuousVectorVariable
		ContinuousVectorVariable* cv = dynamic_cast<ContinuousVectorVariable*>(v);
		if(cv!=0) {
			new ContinuousVectorVariable2SymmetricMatrixVariable(sattr[0],getDAG());
			getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
			return;
		}
		// Try other conversions here, in decreasing order of preference
		string errMsg = "to_symmetric_matrix_XMLParser: could not convert variable " + sattr[1] + " of type " + v->type() + " to desired type";
		error(errMsg.c_str());
	}
	
	// TRANSFORMATIONS: OTHER

	centre_vector_transform_XMLParser::centre_vector_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<centre_vector_transform_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 2;
		const char* attrNames[nattr] = {"id","x"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Obtain random variables
		Variable* vec = getDAG()->get_variable(sattr[1]);
		// Check they are of correct type
		ContinuousVectorVariable* cvec = dynamic_cast<ContinuousVectorVariable*>(vec);
		if(cvec==0) {
			string errMsg = "centre_vector_transform_XMLParser: could not convert variable " + sattr[2] + " of type " + vec->type() + " to continuous vector";
			error(errMsg.c_str());
		}
		new CentreVectorTransform(sattr[0],getDAG());
		getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
		return;
	}
	
	compute_kinship_matrix_transform_XMLParser::compute_kinship_matrix_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<compute_kinship_matrix_transform_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 3;
		const char* attrNames[nattr] = {"id","bip","position"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Obtain random variable
		Variable* v = getDAG()->get_variable(sattr[1]);
		// Try to convert from MatrixVariable
		MatrixVariable* cv = dynamic_cast<MatrixVariable*>(v);
		if(cv!=0) {
			if(sattr[2]=="") {
				new ComputeKinshipMatrixTransform(sattr[0],getDAG());
				getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
				return;
			} else {
				// Convert position to int
				int int_position;
				if(!from_string<int>(int_position,sattr[2])) {
					error("compute_kinship_matrix_transform_XMLParser: could not convert length to int");
				}
				// Create a dummy variable with an internally-generated name. Since it won't be possible to refer to this variable by name, force it to be constant
				// It will not be possible to put a distribution on a constant created this way
				stringstream name;
				name << "_" << sattr[0];
				new DiscreteRV(name.str(),getDAG(),int_position);
				getDAG()->set_constant(name.str());
				new ComputeKinshipMatrixPositionTransform(sattr[0],getDAG());
				getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
				getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],name.str());
				return;
			}
		}
		// Try other conversions here, in decreasing order of preference
		string errMsg = "compute_kinship_matrix_transform_XMLParser: could not convert variable " + sattr[1] + " of type " + v->type() + " to matrix";
		error(errMsg.c_str());
	}
	
	continuous_vector_mixture_component_loglikelihood_XMLParser::continuous_vector_mixture_component_loglikelihood_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<continuous_vector_mixture_component_loglikelihood_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 3;
		const char* attrNames[nattr] = {"id","distribution","rv"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		new ContinuousVectorMixtureComponentLogLikelihood(sattr[2],sattr[1],sattr[0],getDAG());
	}
	
	mean_vector_transform_XMLParser::mean_vector_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<mean_vector_transform_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 2;
		const char* attrNames[nattr] = {"id","x"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Obtain random variables
		Variable* vec = getDAG()->get_variable(sattr[1]);
		// Check they are of correct type
		ContinuousVectorVariable* cvec = dynamic_cast<ContinuousVectorVariable*>(vec);
		if(cvec==0) {
			string errMsg = "mean_vector_transform_XMLParser: could not convert variable " + sattr[2] + " of type " + vec->type() + " to continuous vector";
			error(errMsg.c_str());
		}
		new MeanVector(sattr[0],getDAG());
		getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
		return;
	}
	
	rescale_vector_transform_XMLParser::rescale_vector_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<rescale_vector_transform_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 3;
		const char* attrNames[nattr] = {"id","scalar","vector"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Obtain random variables
		Variable* sca = getDAG()->get_variable(sattr[1]);
		Variable* vec = getDAG()->get_variable(sattr[2]);
		// Check they are of correct type
		ContinuousVariable* csca = dynamic_cast<ContinuousVariable*>(sca);
		ContinuousVectorVariable* cvec = dynamic_cast<ContinuousVectorVariable*>(vec);
		if(csca==0) {
			string errMsg = "rescale_vector_transform_XMLParser: could not convert variable " + sattr[1] + " of type " + sca->type() + " to continuous scalar";
			error(errMsg.c_str());
		}
		if(cvec==0) {
			string errMsg = "rescale_vector_transform_XMLParser: could not convert variable " + sattr[2] + " of type " + vec->type() + " to continuous vector";
			error(errMsg.c_str());
		}
		new RescaleVectorTransform(sattr[0],getDAG());
		getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
		getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
		return;
	}
	
	rescale_symmetric_matrix_transform_XMLParser::rescale_symmetric_matrix_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<rescale_symmetric_matrix_transform_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 3;
		const char* attrNames[nattr] = {"id","scalar","matrix"};
		vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
		// Obtain random variables
		Variable* sca = getDAG()->get_variable(sattr[1]);
		Variable* mat = getDAG()->get_variable(sattr[2]);
		// Check they are of correct type
		ContinuousVariable* csca = dynamic_cast<ContinuousVariable*>(sca);
		SymmetricMatrixVariable* cmat = dynamic_cast<SymmetricMatrixVariable*>(mat);
		if(csca==0) {
			string errMsg = "rescale_symmetric_matrix_transform_XMLParser: could not convert variable " + sattr[1] + " of type " + sca->type() + " to continuous scalar";
			error(errMsg.c_str());
		}
		if(cmat==0) {
			string errMsg = "rescale_symmetric_matrix_transform_XMLParser: could not convert variable " + sattr[2] + " of type " + mat->type() + " to symmetric matrix";
			error(errMsg.c_str());
		}
		new RescaleSymmetricMatrixTransform(sattr[0],getDAG());
		getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
		getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
		return;
	}
	
	symmetric_matrix_sum_transform_XMLParser::symmetric_matrix_sum_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<symmetric_matrix_sum_transform_XMLParser>(master_parser,parent_parser) {
		// Read in the attributes
		const int nattr = 1;
		const char* attrNames[nattr] = {"id"};
		sattr = attributesToStrings(nattr,attrNames,attrs);
		// Don't instantiate the variable until the values have been read in
	}
	
	void symmetric_matrix_sum_transform_XMLParser::implement_characters(const XMLCh* const chars, const XMLSize_t length) {
		string message = "";
		int i;
		for(i=0;i<length;i++) message += chars[i];
		vector<string> operands;
		if(!string_to_vector(operands,message)) error("symmetric_matrix_sum_transform_XMLParser: error interpretting input");
		if(operands.size()==0) error("symmetric_matrix_sum_transform_XMLParser: no operands entered");
		new SymmetricMatrixSumTransform(operands.size(),sattr[0],getDAG());
		for(i=0;i<operands.size();i++) {
			stringstream s;
			s << "operand" << i;
			getDAG()->assign_parameter_to_transformation(sattr[0],s.str(),operands[i]);
		}
	}	
	
	int _GCAT_LMM_LIBRARY_IS_LOADED = 0;
	
	xsd_string load_gcat_lmm_library() {
		if(_GCAT_LMM_LIBRARY_IS_LOADED!=0) {
			throw std::runtime_error("load_gcat_lmm_library(): library already loaded");
		} else {
			_GCAT_LMM_LIBRARY_IS_LOADED = 1;
		}
		// GSL is used, so must set this
		gsl_set_error_handler_off();
		// DISTRIBUTIONS
		distributions_XMLParser::add_child("continuous_vector_mixture",&continuous_vector_mixture_distribution_XMLParser::factory);
		distributions_XMLParser::add_child("marginal_normal_distribution",&marginal_normal_distribution_XMLParser::factory);
		distributions_XMLParser::add_child("multivariate_normal_distribution",&multivariate_normal_distribution_XMLParser::factory);
		// RANDOM VARIABLES: DATA
		data_XMLParser::add_child("continuous_vector_file",&continuous_vector_file_XMLParser::factory);
		data_XMLParser::add_child("continuous_matrix_file",&continuous_matrix_file_XMLParser::factory);
		// RANDOM VARIABLES: PARAMETERS
		parameters_XMLParser::add_child("continuous_vector_file",&continuous_vector_file_XMLParser::factory);		
		// TRANSFORMATIONS: CONVERSIONS
		transformations_XMLParser::add_child("to_diagonal_matrix",&to_diagonal_matrix_XMLParser::factory);
		transformations_XMLParser::add_child("to_symmetric_matrix",&to_symmetric_matrix_XMLParser::factory);
		// TRANSFORMATIONS: OTHER
		transformations_XMLParser::add_child("centre_vector_transform",&centre_vector_transform_XMLParser::factory);
		transformations_XMLParser::add_child("compute_kinship_matrix",&compute_kinship_matrix_transform_XMLParser::factory);
		transformations_XMLParser::add_child("continuous_vector_mixture_component_loglikelihood",&continuous_vector_mixture_component_loglikelihood_XMLParser::factory);
		transformations_XMLParser::add_child("mean",&mean_vector_transform_XMLParser::factory);
		transformations_XMLParser::add_child("rescale_vector_transform",&rescale_vector_transform_XMLParser::factory);
		transformations_XMLParser::add_child("rescale_symmetric_matrix_transform",&rescale_symmetric_matrix_transform_XMLParser::factory);
		transformations_XMLParser::add_child("symmetric_matrix_sum_transform",&symmetric_matrix_sum_transform_XMLParser::factory);
		// SCHEMA
		string s(gcat_lmm1_0_xsd_len,' ');
		unsigned int i;
		for(i=0;i<gcat_lmm1_0_xsd_len;i++) s[i] = gcat_lmm1_0_xsd[i];
		return s;
	}
	
} // namespace gcat_lmm
