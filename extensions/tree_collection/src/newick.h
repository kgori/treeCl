/*
 * newick.h
 *
 *  Created on: May 21, 2010
 *      Author: sadam
 */

#ifndef NEWICK_H_
#define NEWICK_H_

#include "PhyTree.h"
#include <iostream>
#include <string>

namespace newick_parser {
    
    class LexerException: public std::exception {
    private:
        std::string error;
    public:
        LexerException(std::string error) throw () {
            this->error = error;
        }
        
        ~LexerException() throw () {
        }
        
        virtual const char *what() const throw () {
            return error.c_str();
        }
    };
    
    class ParserException: public std::exception {
    private:
        std::string error;
    public:
        ParserException(std::string error) throw () {
            this->error = error;
        }
        
        ~ParserException() throw () {
        }
        
        virtual const char *what() const throw () {
            return error.c_str();
        }
    };
    
    PhyTree::TreePtr parse_newick(std::istream *in) throw (ParserException, LexerException);
    
}
#endif /* NEWICK_H_ */
