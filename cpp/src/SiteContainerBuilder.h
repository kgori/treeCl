/*
 * SiteContainerBuilder.h
 *
 *  Created on: Jan 9, 2014
 *      Author: kgori
 */

#ifndef SITECONTAINERBUILDER_H_
#define SITECONTAINERBUILDER_H_

#include <memory>
#include <stdexcept>
#include <Bpp/Seq/Container.all>
using namespace bpp;
using namespace std;

namespace treeCl {

class SiteContainerBuilder {
public:
	SiteContainerBuilder();
	virtual ~SiteContainerBuilder();
	static shared_ptr<VectorSiteContainer> read_alignment(string filename, string file_format, string datatype, bool interleaved=true) throw (invalid_argument);
private:
	static bool asking_for_fasta(string file_format);
	static bool asking_for_phylip(string file_format);
	static bool asking_for_dna(string datatype);
	static bool asking_for_protein(string datatype);
	static shared_ptr<VectorSiteContainer> read_fasta_dna_file(string filename);
	static shared_ptr<VectorSiteContainer> read_fasta_protein_file(string filename);
	static shared_ptr<VectorSiteContainer> read_phylip_dna_file(string filename, bool interleaved);
	static shared_ptr<VectorSiteContainer> read_phylip_protein_file(string filename, bool interleaved);
};

} /* namespace treeCl */

#endif /* SITECONTAINERBUILDER_H_ */
