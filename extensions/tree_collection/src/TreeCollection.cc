#include <fstream>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include "MinSqTree.h"
#include "ProblemParser.h"

using namespace std;


bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}


int main(int argc, char **argv) {
    int c;

    std::vector<MinSquareTreeCollection::DblMatrix> pmatrices;
    MinSquareTreeCollection::IntMatrix pmapping;
    std::vector<std::string> plabels;
    PhyTree *ptree = NULL;

    MinSquareTreeCollection *mstc = NULL;

    std::string outfilename("out");

    ifstream matrices_data, mapping_data, labels_data, tree_data;
    ofstream statoutfile, treeoutfile;
    int iter(3);

    while ((c=getopt(argc,argv,"D:M:L:T:I:O:")) != EOF) {
        switch (c) {
            case 'D':
                //printf("D = %s\n",optarg);
                if (fileExists(optarg))
                    matrices_data.open(optarg);
                else {
                    cerr << "problem reading " << optarg << endl;
                    abort();
                }
                break;
            case 'M':
                //printf("M = %s\n",optarg);
                if (fileExists(optarg))
                    mapping_data.open(optarg);
                else {
                    cerr << "problem reading " << optarg << endl;
                    abort();
                }
                break;
            case 'L':
                //printf("L = %s\n",optarg);
                if (fileExists(optarg))
                    labels_data.open(optarg);
                else {
                    cerr << "problem reading " << optarg << endl;
                    abort();
                }
                break;
            case 'T':
                //printf("T = %s\n",optarg);
                if (fileExists(optarg))
                    tree_data.open(optarg);
                else {
                    cerr << "problem reading " << optarg << endl;
                    abort();
                }
                break;
            case 'I':
                try {
                    iter = atoi(optarg);
                }
                catch(std::exception &e) {
                    cerr << e.what() << ": couldn't convert " << optarg << " to an integer" << endl;
                    abort();
                }
                if (iter < 1) {
                    cout << "Minimum number of iterations is 1" << endl;
                    iter = 1;
                }
                else if (iter > 100) {
                    cout << "Maximum number of iterations is 100" << endl;
                    iter = 100;
                }
                break;
            case 'O':
                outfilename.assign(optarg);
                break;
            case '?':
            default:
                cerr << "unknown argument:" << optarg << endl;
                exit(1);
                break;
         }
    }

    /* Verify that all files are open  */
    if (!matrices_data.is_open() || !mapping_data.is_open() || !labels_data.is_open()
                               || !tree_data.is_open()) {
       cerr << "usage: TreeCollection -D matrices_file -M mapping_file -L label_file -T initial_topology_file -I number_of_iterations_to_perform -O outfile_basename" << endl;
       exit(1);
    }
    /* Check that everything is ok with output files */
    if (fileExists(outfilename+std::string(".txt"))) {
        cerr << "Error: file" << optarg << ".txt exists" << endl;
        abort();
    }
    if (fileExists(outfilename+std::string(".nwk"))) {
        cerr << "Error: file " << optarg << ".nwk exists" << endl;
        abort();
    }
    treeoutfile.open(outfilename+std::string(".nwk"));
    statoutfile.open(outfilename+std::string(".txt"));
    if (!statoutfile.is_open())  {
        cerr << "problem writing " << outfilename << ".txt. Lack of write permission?" << endl;
        abort();
    }
    else if (!treeoutfile.is_open()) {
        cerr << "problem writing " << outfilename << ".nwk. Lack of write permission?" << endl;
        abort();
    }


    try {
       pmatrices = ProblemParser::parse_matrices(matrices_data);
       pmapping = ProblemParser::parse_mapping(mapping_data);
       plabels = ProblemParser::parse_labels(labels_data);
       ptree = ProblemParser::parse_tree(tree_data);

       mstc = new MinSquareTreeCollection(pmatrices,pmapping,plabels,*ptree);
       delete ptree;
       ptree = NULL;

       cout << "Running " << iter << " iterations" << endl;
       mstc->compute(false, iter);
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

    treeoutfile << mstc->getTree() << endl;
    treeoutfile.close();
    statoutfile << mstc->getScore() << endl;
    statoutfile.close();

    delete mstc;
    return 0;
}
