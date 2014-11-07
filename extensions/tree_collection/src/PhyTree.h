/*
 * Tree.h
 *
 *  Created on: 23 Mar 2010
 *      Author: sadam
 */

#ifndef PHYTREE_H_
#define PHYTREE_H_

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>
#include <sstream>
#include <string>


class PhyTree {
    
public:
    typedef std::shared_ptr<PhyTree> TreePtr;
    typedef std::vector<TreePtr> TreePtrVec;
    typedef PhyTree* RawPtr;
    
private:
    std::vector<TreePtr> children;
    RawPtr parent;
    double branch_length;
    std::string name;
    
    
    void print_prefix(bool affix_comma, bool affix_semicolon, std::string& output) const {
        if (this->isLeaf()) {
            output += this->name + ":" + std::to_string(this->branch_length);
            if (affix_comma)
                output += ",";
        }
        else {
            output += "(";
            for(TreePtrVec::const_iterator i=this->children.begin(); i < this->children.end() - 1; ++i) {
                (*i)->print_prefix(true, false, output);
            }
            (*(this->children.end() - 1))->print_prefix(false, false, output);
            output += "):" + std::to_string(this->branch_length);
            if (affix_comma) output += ",";
        }
        if (affix_semicolon) output += ";";
    }
public:
    void formatNewickM(std::stringstream& newick) const {
        if (this->isLeaf())
            newick << this->getName() << ":" << this->branch_length;
        else {
            newick << "(";
            TreePtrVec::const_iterator i = this->children.begin();
            for(; i < this->children.end() - 1; ++i) {
                (*i)->formatNewickM(newick);
                newick << ",";
            }
            (*i)->formatNewickM(newick);
            newick << "):" << this->branch_length;
        }
    }
    
    std::string formatNewickR() const {
        if(this->isLeaf()) {
            return this->getName();
        } else {
            std::stringstream newick;
            newick << "(";
            TreePtrVec::const_iterator i=this->children.begin();
            newick << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
            for(++i; i < this->children.end(); ++i) {
                newick << "," << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
            }
            newick << ")";
            return newick.str();
        }
    }
    
    
    PhyTree(std::string name="") {
        this->parent = NULL;
        this->branch_length = 0;
        this->name = name;
    }
    
    TreePtr copy() {
        auto out = std::make_shared<PhyTree>();
        out->branch_length = this->branch_length;
        out->name = this->name;
        
        for(TreePtrVec::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            out->addChild((*i)->copy(),(*i)->branch_length);
        }
        return out;
    }
    
    void addChild(TreePtr child, double branch_length = 0) {
        assert(child.get() != this);
        assert(child->parent == nullptr);
        
        this->children.push_back(child);
        child->parent = this;
        child->branch_length = branch_length;
    }
    
    unsigned int indexOf() {
        RawPtr parent = this->parent;
        assert(parent != NULL);
        
        for(unsigned int i=0; i < parent->children.size(); ++i) {
            if(parent->children[i].get() == this) {
                return i;
            }
        }
        
        assert(false);
        return -1;
    }
    
    void pluck() {
        assert(this->parent != NULL);
        
        int index = this->indexOf();
        TreePtrVec::iterator iter = this->parent->children.begin()+index;
        
        this->parent->children.erase(iter);
        this->parent = NULL;
        
        this->branch_length = 0;
    }
    
    TreePtr pluckChild(unsigned int index) {
        TreePtrVec::iterator iter = this->children.begin()+index;
        TreePtr child = *iter;
        
        this->children.erase(iter);
        child->parent = NULL;
        child->branch_length = 0;
        
        return child;
    }
    
    void deleteChild(unsigned int index) {
        TreePtr child = this->pluckChild(index);
        child.reset();
    }
    
    std::string getName() const {
        return this->name;
    }
    
    RawPtr getParent() {
        return this->parent;
    }
    
    const RawPtr getParent() const {
        return this->parent;
    }
    
    double getBranchLength() const {
        return this->branch_length;
    }
    
    PhyTree& operator[](unsigned int i) {
        return *this->children[i];
    }
    
    const PhyTree& operator[](unsigned int i) const {
        return *this->children[i];
    }
    
    TreePtrVec::iterator firstChild() {
        return this->children.begin();
    }
    
    TreePtrVec::iterator lastChild() {
        return this->children.end();
    }
    
    TreePtrVec::const_iterator firstChild() const {
        return this->children.begin();
    }
    
    TreePtrVec::const_iterator lastChild() const {
        return this->children.end();
    }
    
    size_t n_children() const {
        return this->children.size();
    }
    
    bool isLeaf() const {
        return this->children.empty();
    }
    
    void print() const {
        std::string s("");
        this->print_prefix(false, true, s);
    }
    
    void print(std::string& output) {
        this->print_prefix(false, true, output);
    }
    
    std::string formatNewick() const {
        return this->formatNewickR() + ";";
    }
};

#endif /* PHYTREE_H_ */
