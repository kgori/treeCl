/*
 * newick.cpp
 *
 *  Created on: 23 Mar 2010
 *      Author: sadam
 */

#include <iostream>
#include <string>
#include <cctype>
#include <sstream>
#include "newick.h"

using namespace std;

namespace newick_parser {
    
    class TokenBuffer {
        
    private:
        istream *in;
        string last;
        bool valid;
        
        string next_token() throw (LexerException) {
            string token = "";
            while (*this->in) {
                char c = this->in->get();
                if (isspace(c)) {
                    continue;
                } else if (c == ',' || c == ':' || c == '(' || c == ')' || c == ';') {
                    return token + c;
                } else if (isalnum(c) || c == '.' || c == '_') {
                    token += c;
                    while (*this->in && (isalnum(this->in->peek()) || this->in->peek() == '.' || this->in->peek() == '_' || this->in->peek() == '-')) {
                        char c = this->in->get();
                        token += c;
                    }
                    return token;
                } else {
                    throw LexerException("Unexpected character: " + c);
                }
            }
            
            throw LexerException("Unexpected EOF or I/O error");
        }
        
    public:
        TokenBuffer(istream *in) {
            this->in = in;
            valid = false;
        }
        
        string peek() throw (LexerException) {
            if (!valid) {
                last = next_token();
            }
            valid = true;
            return last;
        }
        
        string next() throw (LexerException) {
            if (!valid) {
                last = next_token();
            }
            valid = false;
            return last;
        }
    };
    
    static double parse_double(const string &str) throw () {
        istringstream ss(str.c_str());
        double out = 0;
        
        ss >> out;
        
        return out;
    }
    
    static PhyTree::TreePtr parse_tree(TokenBuffer &buffer) throw (ParserException,
    LexerException) {
        PhyTree::TreePtr t = make_shared<PhyTree>();
        
        string tok = buffer.next();
        
        if (tok != "(")
            throw ParserException("Unexpected token: '" + tok + "', expected: '('");
            
            do {
                PhyTree::TreePtr child = make_shared<PhyTree>();
                tok = buffer.peek();
                
                if (tok == "(")
                    child = parse_tree(buffer);
                    else
                        child = make_shared<PhyTree>(buffer.next());
                        
                        tok = buffer.next();
                        //FIXME
                        if (tok != ":") {
                            double support = parse_double(tok);
                            tok = buffer.next();
                            (void)support;
                        }
                if (tok != ":")
                    throw ParserException("Unexpected token: '" + tok + "', expected: ':'");
                    
                    tok = buffer.next();
                    double dist = parse_double(tok);
                    
                    t->addChild(child, dist);
                    
                    tok = buffer.peek();
                    if (tok == ")") {
                        tok = buffer.next();
                        break;
                    }
                
                tok = buffer.next();
                if (tok != ",")
                    throw ParserException("Unexpected token: '" + tok + "', expected: ','");
                    } while (true);
        
        return t;
    }
    
    PhyTree::TreePtr parse_newick(istream *in) throw (ParserException, LexerException) {
        TokenBuffer buffer(in);
        PhyTree::TreePtr t = parse_tree(buffer);
        
        string tok = buffer.next();
        if (tok == ":") {
            tok = buffer.next();
            (void)parse_double(tok);
            tok = buffer.next();
        }
        
        if (tok != ";")
            throw ParserException("Unexpected token: " + tok);
            
            return t;
    }
    
}
