/*
 * Seq.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: kgori
 */

#include "Distance.h"

#include <sys/_types/_size_t.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "SiteContainerBuilder.h"
#include "ModelFactory.h"

#include "Bpp/Seq/Alphabet.all"
#include "Bpp/Seq/Io/Fasta.h"
#include "Bpp/Seq/Io/Phylip.h"
#include "Bpp/Phyl/Likelihood.all"
#include "Bpp/Numeric/Prob/ConstantDistribution.h"
#include "Bpp/Numeric/Prob/GammaDiscreteDistribution.h"
#include "Bpp/Phyl/Distance/DistanceEstimation.h"
#include "Bpp/Phyl/Distance/BioNJ.h"
#include "Bpp/Phyl/OptimizationTools.h"
#include "Bpp/Phyl/Io/Newick.h"
#include "Bpp/Seq/SymbolListTools.h"
#include "Bpp/Phyl/Simulation.all"

using namespace bpp;
using namespace std;

namespace treeCl {

/* HELPER FUNCTIONS */

template<typename T>
vector<T> merge_vectors(vector<T> &v1, vector<T> &v2) {
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    vector<T> result;
    set_union(v1.begin(), v1.end(), v2.begin(), v2.end(),
              back_inserter(result));
    return result;
}

shared_ptr<VectorSiteContainer> merge(VectorSiteContainer s1, VectorSiteContainer s2)
throw (AlphabetMismatchException, Exception)
{
    if (s1.getAlphabet()->getAlphabetType() != s2.getAlphabet()->getAlphabetType())
        throw AlphabetMismatchException("SiteContainerTools::merge.", s1.getAlphabet(), s2.getAlphabet());

    auto alphabet = s1.getAlphabet();
    const int unknown_char = alphabet->getUnknownCharacterCode();
    size_t l1 = s1.getNumberOfSites();
    size_t l2 = s2.getNumberOfSites();
    auto n1 = s1.getSequencesNames();
    auto n2 = s2.getSequencesNames();
    vector<string> allnames = merge_vectors(n1, n2);

    unique_ptr<VectorSequenceContainer> tmp_container(new VectorSequenceContainer(alphabet));
    for (auto &n : allnames) {
        vector<int> seqvec1;
        vector<int> seqvec2;
        seqvec1.reserve(l1+l2);
        seqvec2.reserve(l2);

        if (s1.hasSequence(n)) {
            auto tmpvec = s1.getSequence(n).getContent();
            seqvec1.insert( seqvec1.end(), tmpvec.begin(), tmpvec.end());
        }
        else {
            vector<int> tmpvec(l1, unknown_char);
            seqvec1.insert( seqvec1.end(), tmpvec.begin(), tmpvec.end());
        }

        if (s2.hasSequence(n)) {
            auto tmpvec = s2.getSequence(n).getContent();
            seqvec2.insert( seqvec2.end(), tmpvec.begin(), tmpvec.end());
        }
        else {
            vector<int> tmpvec(l2, unknown_char);
            seqvec2.insert( seqvec2.end(), tmpvec.begin(), tmpvec.end());
        }
        seqvec1.insert(seqvec1.end(), seqvec2.begin(), seqvec2.end());
        unique_ptr<BasicSequence> seq(new BasicSequence(n, seqvec1, alphabet));

        tmp_container->addSequence(*seq, false);
    }

    shared_ptr<VectorSiteContainer> output(new VectorSiteContainer(*tmp_container));
    return output;
}

Distance::Distance() {
    distances = make_shared<DistanceMatrix>(0);
}

Distance::Distance(string filename, string file_format, string datatype, bool interleaved) {
    read_alignment(filename, file_format, datatype, interleaved);
    fast_compute_distances();
    bionj();
    if (is_dna()) set_model("JCnuc");
    else set_model("WAG01");
    set_gamma_rates();
    set_likelihood_object();
}

Distance::Distance(shared_ptr<VectorSiteContainer>& managed_sequences) {
    shared_ptr<VectorSiteContainer> new_sequences(new VectorSiteContainer(*managed_sequences));
    sequences = new_sequences;
    string type = sequences->getAlphabet()->getAlphabetType();
    if (type == "DNA alphabet") {
        set_dna();
    }
    else if (type == "Proteic alphabet") {
        set_protein();
    }
}

// copy constructor
Distance::Distance(const Distance& other) {
    shared_ptr<VectorSiteContainer> managed_sites(new VectorSiteContainer(*other.sequences));
    sequences = managed_sites;
    string type = sequences->getAlphabet()->getAlphabetType();
    if (type == "DNA alphabet") {
        set_dna();
    }
    else if (type == "Proteic alphabet") {
        set_protein();
    }
    this->model = Distance::copy_model(other.model);
    this->rates = Distance::copy_rates(other.rates);
    if (other.tree_)
        this->tree_ = make_shared<TreeTemplate<Node>>(*other.tree_);
}

// assignment operator
//Distance& Distance::operator=(Distance rhs){
//  sequences = rhs.sequences;
//  tree_ = make_shared<TreeTemplate<Node>>(*rhs.tree_);
//  model = Distance::copy_model(rhs.model);
//  rates = Distance::copy_rates(rhs.rates);
//  distances = make_shared<DistanceMatrix>(*rhs.distances);
//  tl = make_shared<NNIHomogeneousTreeLikelihood>(*rhs.tl);
//  return *this;
//}

Distance Distance::operator+(const Distance& other) {
    shared_ptr<VectorSiteContainer> new_sites = merge(*sequences, *other.sequences);
    return Distance(new_sites);
}

Distance& Distance::operator+=(const Distance& other) {
    shared_ptr<VectorSiteContainer> new_sites = merge(*sequences, *sequences);
    sequences = new_sites;
    return *this;
}

void Distance::read_alignment(string filename, string file_format, string datatype, bool interleaved) {
    sequences = SiteContainerBuilder::read_alignment(filename, file_format, datatype, interleaved);
    string type = sequences->getAlphabet()->getAlphabetType();
    if (type == "DNA alphabet") {
        set_dna();
    }
    else if (type == "Proteic alphabet") {
        set_protein();
    }
    else {
        cout << "Type = " << type << endl;
    }
}

/*Distance Distance::bootstrap_sample() {
    auto new_sequences = SiteContainerTools::bootstrapSites(*sequences.get());
    new_sequences->reindexSites();
    shared_ptr<VectorSiteContainer> managed_new_sequences(new VectorSiteContainer(*new_sequences));
    delete new_sequences;
    Distance ret{managed_new_sequences};
    return ret;
}*/

/*void Distance::write_fasta(string filename) {
    Fasta writer;
    writer.writeAlignment(filename, *sequences);
}*/

/*void Distance::write_phylip(string filename, bool interleaved) {
    Phylip writer{true, !interleaved, 100, true, "  "};
    writer.writeAlignment(filename, *sequences, true);
}*/

void Distance::set_model(string model_name) {
    unique_ptr<ModelFactory> factory(new ModelFactory());
    model = factory->create(model_name);
}

shared_ptr<TreeTemplate<Node>> Distance::string_to_tree(string newick) {
    Newick treeReader;
    stringstream ss;
    ss << newick;
    TreeTemplate<Node> *tmp = treeReader.read(ss);
    auto tree = make_shared<TreeTemplate<Node>>(*tmp);
    delete tmp;
    return tree;
}

void Distance::set_tree(string newick) {
    tree_ = string_to_tree(newick);
    set_likelihood_object();
}

/*void Distance::set_likelihood_object() {
    auto sites_ = sequences->clone();
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    tl = make_shared<NNIHomogeneousTreeLikelihood>(*tree_, *sites_, model.get(), rates.get(), true, false);
    tl->initialize();
    delete sites_;
}*/

string Distance::newick() {
    string ret;
    if (tree_) {
        Newick treeWriter;
        ostringstream s;
        treeWriter.write(*tree_, s);
        ret = s.str();
    }
    return ret;
}

bool Distance::is_dna() {
    return dna && !protein;
}

bool Distance::is_protein() {
    return protein && !dna;
}

void Distance::set_dna() {
    dna = true;
    protein = false;
}

void Distance::compute_distances() {
    size_t n = sequences->getNumberOfSequences();
    vector<string> names = sequences->getSequencesNames();

    for (size_t i = 0; i < n; i++) {
        (*distances)(i, i) = 0;
        for (size_t j = i + 1; j < n; j++) {
            TwoTreeLikelihood* lik =
                new TwoTreeLikelihood(names[i], names[j], *sequences, model.get(), rates.get(), false);
            lik->initialize();
            lik->enableDerivatives(true);
            size_t d = SymbolListTools::getNumberOfDistinctPositions(sequences->getSequence(i), sequences->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sequences->getSequence(i), sequences->getSequence(j));
            lik->setParameterValue("BrLen", g == 0 ? lik->getMinimumBranchLength() : std::max(lik->getMinimumBranchLength(), static_cast<double>(d) / static_cast<double>(g)));
            // Optimization:
            ParameterList params = lik->getBranchLengthsParameters();
            OptimizationTools::optimizeNumericalParameters(lik, params, 0, 1, 0.000001, 1000000, NULL, NULL, false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT);
            // Store results:
            (*distances)(i, j) = lik->getParameterValue("BrLen");
            (*distances)(j, i) = 1.0 / lik->d2f("BrLen", params);
            cout << "i j dist var : " << i << " " << j << " " << (*distances)(i, j) << " " << (*distances)(j, i) << endl;
            delete lik;
        }
    }
}

void Distance::fast_compute_distances() {
    unsigned int s;
    if (is_dna()) {
        s = 4;
    }
    if (is_protein()) {
        s = 20;
    }
    size_t n = sequences->getNumberOfSequences();
    vector<string> names = sequences->getSequencesNames();
    if (distances) distances.reset();
    distances = make_shared<DistanceMatrix>(names);
    for (size_t i = 0; i < n; i++) {
        (*distances)(i, i) = 0;
        for (size_t j=i+1; j < n; j++) {
            size_t d = SymbolListTools::getNumberOfDistinctPositions(sequences->getSequence(i), sequences->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sequences->getSequence(i), sequences->getSequence(j));
            (*distances)(i, j) = jcdist(d, g, s);
            (*distances)(j, i) = jcvar(d, g, s);
        }
    }
}

double Distance::jcdist(size_t d, size_t g, unsigned int s) {
    double p;
    double dist;
    g > 0 ? p = d / g : p = 0;
    dist = - ((s-1)/s) * log(1 - (s/(s-1) * p));
    return dist;
}

double Distance::jcvar(size_t d, size_t g, unsigned int s) {
    double p;
    double var;
    double denom;
    g > 0 ? p = d / g : p = 0;
    denom = 1 - ((s / (s-1)) * p);
    g > 0 ? var = (p * (1 - p)) / (g * denom * denom): var = 0;
    return var
}

/*void Distance::optimize_numerical(double tolerance, long max_calls) {
    ParameterList tl_params = tl->getParameters();
    OptimizationTools::optimizeNumericalParameters2(tl.get(), tl_params, 0, tolerance, max_calls, NULL, NULL, false, false, 0);
    TreeTemplate<Node> optimized_tree = tl->getTree();
    tree_ = make_shared<TreeTemplate<Node>>(optimized_tree);
}*/

/*void Distance::optimize_tree(double tolerance_before, double tolerance_during, long max_calls, int nnis_per_round) {
    ParameterList tl_params = tl->getParameters();
    auto new_tl = OptimizationTools::optimizeTreeNNI2(tl.get(), tl_params, true, tolerance_before, tolerance_during, max_calls, nnis_per_round, NULL, NULL, false, 0);
    tl = make_shared<NNIHomogeneousTreeLikelihood>(*new_tl);
//  delete new_tl;
    TreeTemplate<Node> optimized_tree = tl->getTree();
    tree_ = make_shared<TreeTemplate<Node>>(optimized_tree);
}*/

Distance Distance::simulate() {
    unique_ptr<HomogeneousSequenceSimulator> simulator(new HomogeneousSequenceSimulator(model.get(), rates.get(), tree_.get()));
    size_t number_of_sites = this->get_alignment_length();
    SiteContainer *tmp = simulator->simulate(number_of_sites);
    shared_ptr<VectorSiteContainer> sim_sites(new VectorSiteContainer(*tmp));
    delete tmp;
    auto ret = Distance(*this);
    auto alphabet = sequences->getAlphabet();
    int gapCode = alphabet->getGapCharacterCode();
    size_t numseq = this->get_number_of_sequences();
    for (unsigned int i = 0; i < numseq; i++) {
        BasicSequence oldseq = sequences->getSequence(i);
        BasicSequence newseq = sim_sites->getSequence(i);
        for (unsigned int j = 0; j < number_of_sites; j++) {
            if (alphabet->isGap(oldseq[j])) {
                newseq.setElement(j, gapCode);
            }
        }
        string name = newseq.getName();
        sim_sites->setSequence(name, newseq, true);
    }
    ret.sequences = sim_sites;
    return ret;
}

void Distance::set_protein() {
    dna = false;
    protein = true;
}

size_t Distance::get_alignment_length() {
    return sequences->getNumberOfSites();
}

size_t Distance::get_number_of_sequences() {
    return sequences->getNumberOfSequences();
}

void Distance::set_gamma_rates(int ncat, double alpha) {
    rates = make_shared<GammaDiscreteDistribution>(ncat, alpha, alpha);
    rates->aliasParameters("alpha", "beta");
}

void Distance::set_constant_rates(double value) {
    rates = make_shared<ConstantDistribution>(value);
}

/*string Distance::bionj(bool fast) {
    if (!distances && fast) fast_compute_distances();
    else if (!distances && !fast) compute_distances();
    BioNJ bionj(*distances, false, true, false); // rooted=false, positiveLengths=true, verbose=false
    TreeTemplate<Node> *njtree = bionj.getTree();
    tree_ = make_shared<TreeTemplate<Node>>(*njtree);
    delete njtree;
    return newick();
}
*/
/*double Distance::get_likelihood() {
    if (!tl) set_likelihood_object();
    return -tl->getValue();
}*/

/*double Distance::get_likelihood(string newick) {
    set_tree(newick);
    double logL = tl->getValue();
    return -logL;
}*/

shared_ptr<SubstitutionModel> Distance::copy_model(shared_ptr<SubstitutionModel> other_model) {
    shared_ptr<SubstitutionModel> model;
    if (other_model) {
        string modelname{other_model->getName()};
        ParameterList modelparams = other_model->getParameters();
        auto factory = new ModelFactory();
        model = factory->create(modelname);
        delete factory;
        model->setAllParametersValues(modelparams);
    }
    return model;
}

shared_ptr<DiscreteDistribution> Distance::copy_rates(
        shared_ptr<DiscreteDistribution> other_rates) {
    shared_ptr<DiscreteDistribution> rates;
    if (other_rates) {
        string ratesname{other_rates->getName()};
        if (ratesname == "Gamma") {
            size_t ncat = other_rates->getNumberOfCategories();
            double alpha = other_rates->getParameterValue("alpha");
            rates = make_shared<GammaDiscreteDistribution>(ncat, alpha, alpha);
            rates->aliasParameters("alpha", "beta");
        }
        else {
            double value = other_rates->getParameterValue("value");
            rates = make_shared<ConstantDistribution>(value);
        }
    }
    return rates;
}

} /* namespace treeCl */
