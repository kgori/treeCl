#include <utility>
#include "../src/MinSqTree.h"
#include "../src/ProblemParser.h"

using namespace std;

double fit(std::string matrices, std::string mapping, std::string labels, std::string tree)
{
    std::vector<MinSquareTreeCollection::DblMatrix> pmatrices;
    MinSquareTreeCollection::IntMatrix pmapping;
    std::vector<std::string> plabels;
    PhyTree::TreePtr ptree;
    MinSquareTreeCollection *mstc = NULL;

    try {
        pmatrices = ProblemParser::parse_matrices(matrices);
        pmapping = ProblemParser::parse_mapping(mapping);
        plabels = ProblemParser::parse_labels(labels);
        ptree = ProblemParser::parse_tree(tree);

        mstc = new MinSquareTreeCollection(pmatrices, pmapping, plabels, *ptree);
    }

    catch(std::exception &e) {
        if(mstc) {
            delete mstc;
        }
        cerr << e.what() << endl;
    }

    double result = mstc->getScore();

    delete mstc;
    return result;
}

std::pair<std::string, double> compute(std::string matrices, std::string mapping, std::string labels, std::string tree,
                                       int iter, bool loglik, bool keep_topology, bool quiet)
{
    std::vector<MinSquareTreeCollection::DblMatrix> pmatrices;
    MinSquareTreeCollection::IntMatrix pmapping;
    std::vector<std::string> plabels;
    PhyTree::TreePtr ptree;

    std::pair<std::string, double> result;

    pmatrices = ProblemParser::parse_matrices(matrices);
    pmapping = ProblemParser::parse_mapping(mapping);
    plabels = ProblemParser::parse_labels(labels);
    ptree = ProblemParser::parse_tree(tree);

    MinSquareTreeCollection mstc(pmatrices, pmapping, plabels, *ptree);
    mstc.compute(keep_topology, iter, quiet);
    mstc.getTree();
    if (loglik)
        result = std::make_pair(mstc.newick, mstc.getLogLikelihood());
    else
        result = std::make_pair(mstc.newick, mstc.getScore());
    return result;
}
