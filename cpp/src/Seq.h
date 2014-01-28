/*
 * Seq.h
 *
 *  Created on: Jan 8, 2014
 *      Author: kgori
 */

#ifndef SEQ_H_
#define SEQ_H_

#include <iostream>
#include <memory>
#include <Bpp/Seq/Container.all>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include "Bpp/Phyl/Likelihood.all"
#include "Bpp/Phyl/BipartitionList.h"
#include "Bpp/Phyl/Node.h"
#include <Bpp/Numeric/Prob.all>
#include <Bpp/Phyl/Tree.h>
#include "Bpp/Seq/Alphabet.all"
#include "Bpp/Seq/Io/Fasta.h"
#include "Bpp/Seq/Io/Phylip.h"
#include "Bpp/Numeric/Prob/ConstantDistribution.h"
#include "Bpp/Numeric/Prob/GammaDiscreteDistribution.h"
#include "Bpp/Phyl/Distance/DistanceEstimation.h"
#include "Bpp/Phyl/Distance/BioNJ.h"
#include "Bpp/Phyl/OptimizationTools.h"
#include "Bpp/Phyl/Io/Newick.h"
#include "Bpp/Seq/SymbolListTools.h"
using namespace std;
using namespace bpp;

namespace treeCl {

    class Seq {
        public :
    		Seq();
    	    Seq(string filename, string file_format, string datatype, bool interleaved=true);
            Seq(shared_ptr<VectorSiteContainer>& managed_sequences);
            Seq(const Seq& other);
            // Seq& operator=(Seq rhs);
            Seq operator+(const Seq& other);
            Seq& operator+=(const Seq& other);
            void read_alignment(string filename, string file_format, string datatype, bool interleaved=true);
            void write_fasta(string filename);
            void write_phylip(string filename, bool interleaved=true);
            size_t get_alignment_length();
            size_t get_number_of_sequences();
            Seq bootstrap_sample();
            void set_model(string model_name);
            void set_gamma_rates(int ncat=4, double alpha=1.0);
            void compute_distances();
            void fast_compute_distances();
            void set_constant_rates(double value=1.0);
            bool is_dna();
            bool is_protein();
            string bionj(bool fast=true);
            void optimize_numerical(double tolerance=10, long max_calls=10000);
            void optimize_tree(double tolerance_before=10, double tolerance_during=100, long max_calls=10000, int nnis_per_round=4);
            void set_tree(string newick);
            double get_likelihood();
            double get_likelihood(string newick);
            void set_likelihood_object();
            Seq simulate();
            string newick();
            shared_ptr<TreeTemplate<Node>> string_to_tree(string newick);
            shared_ptr<VectorSiteContainer> sequences;
            shared_ptr<SubstitutionModel> model;
            shared_ptr<DiscreteDistribution> rates;
            shared_ptr<DistanceMatrix> distances;
            shared_ptr<Tree> tree_;
            shared_ptr<NNIHomogeneousTreeLikelihood> tl;
        private :
            static shared_ptr<SubstitutionModel> copy_model(shared_ptr<SubstitutionModel> other_model);
            static shared_ptr<DiscreteDistribution> copy_rates(shared_ptr<DiscreteDistribution> other_rates);
            double jcdist(double d, double g, double s);
            bool dna{false};
            bool protein{false};
            void set_dna();
            void set_protein();

    };
} /* namespace treeCl */

#endif /* SEQ_H_ */
