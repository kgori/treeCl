#include <utility>
#include "../src/MinSqTree.h"
#include "../src/ProblemParser.h"

using namespace std;

double fit(std::string matrices, std::string mapping, std::string labels, std::string tree)
{
    std::vector<MinSquareTreeCollection::DblMatrix> pmatrices;
    MinSquareTreeCollection::IntMatrix pmapping;
    std::vector<std::string> plabels;
    PhyTree *ptree = NULL;
    MinSquareTreeCollection *mstc = NULL;

    try {
        pmatrices = ProblemParser::parse_matrices(matrices);
        pmapping = ProblemParser::parse_mapping(mapping);
        plabels = ProblemParser::parse_labels(labels);
        ptree = ProblemParser::parse_tree(tree);

        mstc = new MinSquareTreeCollection(pmatrices, pmapping, plabels, *ptree);
        delete ptree;
        ptree = NULL;
    }

    catch(std::exception &e) {
        if(ptree) {
            delete ptree;
        }
        if(mstc) {
            delete mstc;
        }
        cerr << e.what() << endl;
    }

    double result = mstc->getScore();

    delete mstc;
    return result;
}

std::pair<std::string, double> compute(std::string matrices, std::string mapping, std::string labels, std::string tree, int iter, bool keep_topology, bool quiet) throw()
{
    std::vector<MinSquareTreeCollection::DblMatrix> pmatrices;
    MinSquareTreeCollection::IntMatrix pmapping;
    std::vector<std::string> plabels;
    PhyTree *ptree = NULL;
    MinSquareTreeCollection *mstc = NULL;

    try {
        pmatrices = ProblemParser::parse_matrices(matrices);
        pmapping = ProblemParser::parse_mapping(mapping);
        plabels = ProblemParser::parse_labels(labels);
        ptree = ProblemParser::parse_tree(tree);

        mstc = new MinSquareTreeCollection(pmatrices, pmapping, plabels, *ptree);
        delete ptree;
        ptree = NULL;

        mstc->compute(keep_topology, iter, quiet);
    }

    catch(std::exception &e) {
        if(ptree) {
            delete ptree;
        }
        if(mstc) {
            delete mstc;
        }
        cerr << e.what() << endl;
    }

    std::string out_tree(mstc->getTree());
    double out_score(mstc->getScore());
    std::pair<std::string, double> result = std::make_pair(out_tree, out_score);

    delete mstc;
    return result;
}
