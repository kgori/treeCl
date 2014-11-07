#include <iostream>
#include <sstream>
#include <string>
#include "Eigen/Core"
#include <cstdio>
#include <vector>

#include "ProblemParser.h"
#include "MinSqTree.h"
#include "newick.h"

using namespace std;

namespace ProblemParser {
    vector<MinSquareTreeCollection::DblMatrix> parse_matrices(istream &matrices) throw (ParseException) {
        vector<MinSquareTreeCollection::DblMatrix> pmatrices;
        int n_matrices=0;
        std::stringstream errMsg;
        
        matrices >> n_matrices;
        for(int k=0; k<n_matrices; ++k) {
            int rows=0, cols=0, index=0;
            matrices >> rows >> cols >> index;
            if(index != k+1) {
                throw ParseException("wrong index");
            }
            MinSquareTreeCollection::DblMatrix m(rows,cols);
            for(int i=0; i<rows; ++i) {
                for(int j=0; j<cols; ++j) {
                    matrices >> m(i,j);
                }
            }
            
            if(!matrices) {
                errMsg << "I/O error parsing matrices (k=" << k << ")";
                throw ParseException(errMsg.str());
                }
                
                pmatrices.push_back(m);
                }
                
                return pmatrices;
                }
                
                MinSquareTreeCollection::IntMatrix parse_mapping(istream &mapping) throw (ParseException) {
                    int n_matrices=0;
                    int n_genomes=0;
                    
                    mapping >> n_matrices;
                    mapping >> n_genomes;
                    MinSquareTreeCollection::IntMatrix pmapping(n_matrices,n_genomes);
                    for(int i=0; i<n_matrices; ++i) {
                        for(int j=0; j<n_genomes; ++j) {
                            mapping >> pmapping(i,j);
                        }
                    }
                    
                    if(!mapping) {
                        throw ParseException("I/O error parsing mapping");
                    }
                    
                    return pmapping;
                }
                
                vector<string> parse_labels(istream &labels) throw (ParseException) {
                    int n_genomes=0;
                    
                    labels >> n_genomes;
                    vector<string> plabels;
                    for(int i=0; i<n_genomes; ++i) {
                        string label;
                        labels >> label;
                        plabels.push_back(label);
                    }
                    
                    if(!labels) {
                        throw ParseException("I/O error parsing labels");
                    }
                    return plabels;
                }
                
                PhyTree::TreePtr parse_tree(istream &tree) throw (ParseException) {
                    PhyTree::TreePtr ptree = newick_parser::parse_newick(&tree);
                    return ptree;
                }
                
                void parse_and_run(istream &matrices, istream &mapping, istream &labels, istream &tree) throw (ParseException) {
                    vector<MinSquareTreeCollection::DblMatrix> pmatrices = parse_matrices(matrices);
                    MinSquareTreeCollection::IntMatrix pmapping = parse_mapping(mapping);
                    vector<string> plabels = parse_labels(labels);
                    PhyTree::TreePtr ptree = parse_tree(tree);
                    
                    MinSquareTreeCollection matc(pmatrices, pmapping, plabels, *ptree);
                    matc.compute(false);
                    
                }
                
                void parse_and_run(const string &matrices, const string &mapping, const string &labels, const string &tree) throw (ParseException) {
                    stringstream smatrices(matrices);
                    stringstream smapping(mapping);
                    stringstream slabels(labels);
                    stringstream stree(tree);
                    parse_and_run(smatrices,smapping,slabels,stree);
                }
                
                vector<MinSquareTreeCollection::DblMatrix> parse_matrices(const string &matrices) throw (ParseException) {
                    stringstream smatrices(matrices);
                    return parse_matrices(smatrices);
                }
                
                MinSquareTreeCollection::IntMatrix parse_mapping(const string &mapping) throw (ParseException) {
                    stringstream smapping(mapping);
                    return parse_mapping(smapping);
                }
                
                vector<string> parse_labels(const string &labels) throw (ParseException) {
                    stringstream slabels(labels);
                    return parse_labels(slabels);
                }
                
                PhyTree::TreePtr parse_tree(const string &tree) throw (ParseException) {
                    stringstream stree(tree);
                    return parse_tree(stree);
                }
                
                }
