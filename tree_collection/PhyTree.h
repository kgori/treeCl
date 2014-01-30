/*
 * Tree.h
 *
 *  Created on: 23 Mar 2010
 *      Author: sadam
 */

#ifndef PHYTREE_H_
#define PHYTREE_H_

#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <sstream>

class PhyTree {
private:
	std::vector<PhyTree*> children;
	PhyTree *parent;
	double branch_length;
	std::string name;

	void print_prefix(std::string prefix) const {
		std::cout << prefix << "(" << this->branch_length << ") " << this->name << std::endl;
		for(std::vector<PhyTree*>::const_iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->print_prefix(prefix+"  ");
		}
	}

	std::string formatNewickR() const {
		if(this->isLeaf()) {
			return this->getName();
		} else {
			std::stringstream newick;
			newick << "(";
			std::vector<PhyTree*>::const_iterator i=this->children.begin();
			newick << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
			for(++i; i < this->children.end(); ++i) {
				newick << "," << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
			}
			newick << ")";
			return newick.str();
		}
	}

public:
	PhyTree(std::string name="") {
		this->parent = NULL;
		this->branch_length = 0;
		this->name = name;
	}

	~PhyTree() {
		assert(this->parent == NULL);
		for(std::vector<PhyTree*>::reverse_iterator i=this->children.rbegin(); i < this->children.rend(); ++i) {
			PhyTree *child = *i;
			child->parent = NULL;
			child->branch_length = 0;

			delete child;
		}
	}

	PhyTree* copy() {
		PhyTree* out = new PhyTree();
		out->branch_length = this->branch_length;
		out->name = this->name;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			out->addChild((*i)->copy(),(*i)->branch_length);
		}
		return out;
	}

	void addChild(PhyTree *child, double branch_length = 0) {
		assert(child != this);
		assert(child->parent == NULL);

		this->children.push_back(child);
		child->parent = this;
		child->branch_length = branch_length;
	}

	unsigned int indexOf() {
		PhyTree *parent = this->parent;
		assert(parent != NULL);

		for(unsigned int i=0; i < parent->children.size(); ++i) {
			if(parent->children[i] == this) {
				return i;
			}
		}

		assert(false);
		return -1;
	}

	void pluck() {
		assert(this->parent != NULL);

		int index = this->indexOf();
		std::vector<PhyTree*>::iterator iter = this->parent->children.begin()+index;

		this->parent->children.erase(iter);
		this->parent = NULL;

		this->branch_length = 0;
	}

	PhyTree* pluckChild(unsigned int index) {
		std::vector<PhyTree*>::iterator iter = this->children.begin()+index;
		PhyTree *child = *iter;

		this->children.erase(iter);
		child->parent = NULL;
		child->branch_length = 0;

		return child;
	}

	void deleteChild(unsigned int index) {
		PhyTree *child = this->pluckChild(index);
		delete child;
	}

	std::string getName() const {
		return this->name;
	}

	PhyTree *getParent() {
		return this->parent;
	}

	const PhyTree *getParent() const {
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

	std::vector<PhyTree*>::iterator firstChild() {
		return this->children.begin();
	}

	std::vector<PhyTree*>::iterator lastChild() {
		return this->children.end();
	}

	std::vector<PhyTree*>::const_iterator firstChild() const {
		return this->children.begin();
	}

	std::vector<PhyTree*>::const_iterator lastChild() const {
		return this->children.end();
	}

	unsigned int n_children() const {
		return this->children.size();
	}

	bool isLeaf() const {
		return this->children.empty();
	}

	void print() const {
		this->print_prefix("");
	}

	std::string formatNewick() const {
		return this->formatNewickR() + ";";
	}
};

#endif /* PHYTREE_H_ */
