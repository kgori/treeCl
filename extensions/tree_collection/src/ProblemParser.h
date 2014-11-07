#include <iostream>
#include <string>
#include <vector>

#include "PhyTree.h"
#include "MinSqTree.h"

namespace ProblemParser {
    class ParseException: public std::exception {
    private:
        std::string error;
    public:
        ParseException(const std::string &error) throw () {
            this->error = "Error while parsing problem: " + error;
        }
        
        ~ParseException() throw () {
        }
        
        virtual const char *what() const throw () {
            return error.c_str();
        }
    };
    
    std::vector<MinSquareTreeCollection::DblMatrix> parse_matrices(const std::string &matrices) throw (ParseException);
    MinSquareTreeCollection::IntMatrix parse_mapping(const std::string &mapping) throw (ParseException);
    std::vector<std::string> parse_labels(const std::string &labels) throw (ParseException);
    PhyTree::TreePtr parse_tree(const std::string &tree) throw (ParseException);
    
    std::vector<MinSquareTreeCollection::DblMatrix> parse_matrices(std::istream &matrices) throw (ParseException);
    MinSquareTreeCollection::IntMatrix parse_mapping(std::istream &mapping) throw (ParseException);
    std::vector<std::string> parse_labels(std::istream &labels) throw (ParseException);
    PhyTree::TreePtr parse_tree(std::istream &tree) throw (ParseException);
    
    void parse_and_run(const std::string &matrices, const std::string &mapping, const std::string lables, const std::string &tree) throw (ParseException);
    void parse_and_run(std::istream &matrices, std::istream &mapping, std::istream &labels, std::istream &tree) throw (ParseException);
}
