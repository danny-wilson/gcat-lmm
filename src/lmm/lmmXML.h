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
#ifndef _GCAT_LMM_XML_H_
#define _GCAT_LMM_XML_H_
#include <DAG/DAGXMLParser.h>
#include <Inference/MCMC/MCMC.h>

using namespace gcat;

namespace gcat_lmm {

	// DISTRIBUTIONS
	
/*	<xs:element name="continuous_vector_mixture" substitutionGroup="abstract_distribution">
		<xs:complexType>
			<xs:simpleContent>
				<xs:extension base="xs:string">
					<xs:attribute name="id" type="xs:string" use="required"/>
					<xs:attribute name="length" type="xs:string" use="required"/>
					<xs:attribute name="p" type="xs:string" default=""/>
				</xs:extension>
			</xs:simpleContent>
		</xs:complexType>
	 </xs:element>
 */	 
	class continuous_vector_mixture_distribution_XMLParser : public DAGXMLParserTemplate<continuous_vector_mixture_distribution_XMLParser> {
		vector<string> sattr;
	public:
		continuous_vector_mixture_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
		void implement_characters(const XMLCh* const chars, const XMLSize_t length);
	};
	
/*	<xs:element name="marginal_normal_distribution" substitutionGroup="abstract_distribution">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="length" type="xs:string" use="required"/>
			<xs:attribute name="cov" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/
	class marginal_normal_distribution_XMLParser : public DAGXMLParserTemplate<marginal_normal_distribution_XMLParser> {
	public:
		marginal_normal_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="multivariate_normal_distribution" substitutionGroup="abstract_distribution">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="length" type="xs:string" use="required"/>
			<xs:attribute name="mean" type="xs:string" use="required"/>
			<xs:attribute name="cov" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/
	class multivariate_normal_distribution_XMLParser : public DAGXMLParserTemplate<multivariate_normal_distribution_XMLParser> {
	public:
		multivariate_normal_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
	// RANDOM VARIABLES
	
/*	<xs:element name="continuous_vector_file" substitutionGroup="abstract_parameter">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" default=""/>
			<xs:attribute name="file" type="xs:string" use="required"/>
			<xs:attribute name="sep" type="xs:string" default=""/>
		</xs:complexType>
	</xs:element>
*/
	class continuous_vector_file_XMLParser : public DAGXMLParserTemplate<continuous_vector_file_XMLParser> {
	public:
		continuous_vector_file_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser, const bool constant=false);
		// Specialized static member functions
		static DAGXMLParser* factory_constant(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) {
			return (DAGXMLParser*)(new continuous_vector_file_XMLParser(uri,localname,qname,attrs,master_parser,parent_parser,true));
		}
	};
	
/*	<xs:element name="continuous_matrix_file" substitutionGroup="abstract_parameter">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" default=""/>
			<xs:attribute name="file" type="xs:string" use="required"/>
			<xs:attribute name="sep" type="xs:string" default=""/>
		</xs:complexType>
	</xs:element>
*/
	class continuous_matrix_file_XMLParser : public DAGXMLParserTemplate<continuous_matrix_file_XMLParser> {
	public:
		continuous_matrix_file_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser, const bool constant=false);
		// Specialized static member functions
		static DAGXMLParser* factory_constant(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) {
			return (DAGXMLParser*)(new continuous_matrix_file_XMLParser(uri,localname,qname,attrs,master_parser,parent_parser,true));
		}
	};

	// TRANSFORMATIONS: CONVERSIONS
	
/*	<xs:element name="to_diagonal_matrix" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="x" type="xs:string" use="required"/>
			<xs:attribute name="length" type="xs:string" default=""/>
		</xs:complexType>
	</xs:element>
*/
	class to_diagonal_matrix_XMLParser : public DAGXMLParserTemplate<to_diagonal_matrix_XMLParser> {
	public:
		to_diagonal_matrix_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="to_symmetric_matrix" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="x" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/	 
	class to_symmetric_matrix_XMLParser : public DAGXMLParserTemplate<to_symmetric_matrix_XMLParser> {
	public:
		to_symmetric_matrix_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
	// TRANSFORMATIONS: OTHER
	
/*	<xs:element name="centre_vector_transform" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="x" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
	 */	 
	class centre_vector_transform_XMLParser : public DAGXMLParserTemplate<centre_vector_transform_XMLParser> {
	public:
		centre_vector_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
		
/*	<xs:element name="compute_kinship_matrix" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="bip" type="xs:string" use="required"/>
			<xs:attribute name="position" type="xs:string" default=""/>
		</xs:complexType>
	 </xs:element>
*/
	class compute_kinship_matrix_transform_XMLParser : public DAGXMLParserTemplate<compute_kinship_matrix_transform_XMLParser> {
	public:
		compute_kinship_matrix_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="mean" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="vector" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/	 
	class mean_vector_transform_XMLParser : public DAGXMLParserTemplate<mean_vector_transform_XMLParser> {
	public:
		mean_vector_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="continuous_vector_mixture_component_loglikelihood">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" use="required"/>
			<xs:attribute name="rv" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
	class continuous_vector_mixture_component_loglikelihood_XMLParser : public DAGXMLParserTemplate<continuous_vector_mixture_component_loglikelihood_XMLParser> {
	public:
		continuous_vector_mixture_component_loglikelihood_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="rescale_vector_transform" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="scalar" type="xs:string" use="required"/>
			<xs:attribute name="vector" type="xs:string" use="required"/>
		</xs:complexType>
	 </xs:element>
*/	 
	class rescale_vector_transform_XMLParser : public DAGXMLParserTemplate<rescale_vector_transform_XMLParser> {
	public:
		rescale_vector_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="rescale_symmetric_matrix_transform" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="scalar" type="xs:string" use="required"/>
			<xs:attribute name="matrix" type="xs:string" use="required"/>
		</xs:complexType>
	 </xs:element>
*/	 
	class rescale_symmetric_matrix_transform_XMLParser : public DAGXMLParserTemplate<rescale_symmetric_matrix_transform_XMLParser> {
	public:
		rescale_symmetric_matrix_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	};
	
/*	<xs:element name="symmetric_matrix_sum_transform" substitutionGroup="abstract_transformation">
		<xs:complexType>
			<xs:simpleContent>
				<xs:extension base="xs:string">
					<xs:attribute name="id" type="xs:string" use="required"/>
				</xs:extension>
			</xs:simpleContent>
		</xs:complexType>
	</xs:element>
*/
	class symmetric_matrix_sum_transform_XMLParser : public DAGXMLParserTemplate<symmetric_matrix_sum_transform_XMLParser> {
		vector<string> sattr;
	public:
		symmetric_matrix_sum_transform_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
		void implement_characters(const XMLCh* const chars, const XMLSize_t length);
	};
	
	
	// Load the library
	xsd_string load_gcat_lmm_library();
	
} // namespace gcat_lmm

#endif//_GCAT_LMM_XML_H_
