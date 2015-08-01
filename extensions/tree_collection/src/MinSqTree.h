#ifndef _MINSQTREE_H
#define _MINSQTREE_H
#include "Eigen/Core"
#include <string>
#include <vector>
#include <exception>

#include "PhyTree.h"

class MinSquareTreeCollection
{
public:
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> DblMatrix;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> IntMatrix;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> DblVector;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> IntVector;
    
    class ParameterException: public std::exception {
    private:
        std::string error;
    public:
        ParameterException(const std::string &error) throw () {
            this->error = "Wrong input parameters: " + error;
        }
        
        ~ParameterException() throw () {
        }
        
        virtual const char *what() const throw () {
            return error.c_str();
        }
    };
    
    class RuntimeException: public std::exception {
    private:
        std::string error;
    public:
        RuntimeException(const std::string &error) throw () {
            this->error = "Runtime error in MinSqTree: " + error;
        }
        
        ~RuntimeException() throw () {
        }
        
        virtual const char *what() const throw () {
            return error.c_str();
        }
    };
    
    
private:
    /** edge description */
    struct edge_t {
        double len;
        double support;
        int From, To;
    };
    
    /** edge description */
    struct edgec_t {
        DblVector len;
        int From, To;
        double len0;
        
        void alloc(int n) {
            len = DblVector::Zero(n);
        }
    };
    
    /* global variables */
    std::vector<DblMatrix> aDistVar;
    IntMatrix aMap;
    int NT;
    double globmin;
    DblVector globminA, tmpA;
    std::vector<edgec_t> EdgeC;
    double MST_Qual = 0;
    int NewInternalNode, NewEdgeIndex;
    
    /** Minimum edge length */
    double MinLen;
    
    /** number of external nodes */
    int ne;
    
    std::vector<edge_t> Edge, BestEdge;
    
    /* external nodes (sequences) are numbered from 0 to ne-1 */
    /* internal nodes (common ancestors) are numbered from ne to 2ne-3 */
    /* edges are numbered from 0 to 2ne-4 */
    
    /* we need to allocate one val per vertex */
    DblMatrix ConShortestPathC;
    IntVector ShortestLabel;
    std::vector<std::string> Labels;
    
    /* internal node edge incidences */
    IntMatrix inc;
    
    /* status */
    bool initialized;
    bool computed;
    
    /* Methods */
    void SolveRestricted( const DblMatrix &Z, const DblVector &B, DblVector &L, int n, int restr );
    int BestRestrictionRemoval( const DblMatrix &Z, const DblVector &B, DblVector &L, int n, int restr );
    double FourSubtree( double Wab, double Mab, double Wac, double Mac, double Wad, double Mad, double Wbc, double Mbc, double Wbd, double Mbd, double Wcd, double Mcd, DblVector &L);
    void printTreeC();
    void printTree();
    void MS_ShortestPathCollection( int from, int ExcludedEdge, int label );
    void MS_ShortestPathOne( int from, int ExcludedEdge, int k, int label );
    double DistanceFitCollection();
    double LogLikelihoodFitCollection();
    void IncidencesC();
    void getSons(int e0, int n, int* e1, int* e2);
    void LabelNonExistEdgesR(int e0, int n_papa);
    void Fix001Case(int er, int e0, int n_papa);
    void LabelNonExistEdges();
    void FindQuartet(int i, int e0, int n_papa, int* at, int* et, int* x, int* ia, int* ix, int* a2, int* ee2, int* ia2, int* ex2);
    void FindSplit(int i, int e0, int n_papa, int* x, int* ix, int* at, int* ia, int *et, int* ex);
    int CountOrLabelPath(int i, int x2, int x1, int ex2, double dL, int* ise0);
    void ThreeOptimSubset(double wab, double mab, double wac, double mac, double wbc, double mbc,double *T);
    int CountOrLabelPathTriplet(int i, int A, int x0, int eA, double dL);
    void FitTriplet(int k, int A, int B, int C, int  eA, int eB, int eC, int x0);
    void FitQuartet(int k, int A, int B, int C, int D, int eA, int eB, int eC, int eD, int x1, int x2, int ex2, int qOk, int e0, int allEdges);
    void GotoLeaf(int i, int e0, int n_papa, int* eL, int *L) ;
    void FitEdge(int i, int e0, int allEdges);
    void FitLabeledEdgesC(int allEdges);
    int FourOptimCollection(int DoNotModify, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W);
    double FiveSubtreeCollection( double Wab, double Mab, double Wac, double Mac, double Wad, double Mad, double Wae, double Mae, double Wbc, double Mbc, double Wbd, double Mbd, double Wbe, double Mbe, double Wcd, double Mcd, double Wce, double Mce, double Wde, double Mde, DblVector &L, int Skip) ;
    void SwapFiveSubtree(int e1, int e2, int e4, int e6, int e7, int ix, int iw, int iy, int e3, int e5);
    int DecodeFlags(int f, int *flag);
    int FiveOptimCollectionExtended345(int e3,int e4,int e5,int vc,int ix, int iy,int iw, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W, int *ptrnswap);
    int ModifySkip(int flag,int p0,int p1,int p2,int p3,int p4);
    void delPathLength(int k, int from, int ExcludedEdge);
    int FiveOptimCollection345(int DoNotModify, int e3,int e4,int e5,int vc,int ix, int iy,int iw, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W, int *ptrnswap);
    int FiveOptimCollection(int mode, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W);
    int MapTree2InternalC( const PhyTree &t ) throw (ParameterException);
    double EdgeCLength( int i );
    PhyTree::TreePtr MST_TreeCR( int node, int edgefrom);
    
public:
    MinSquareTreeCollection( const std::vector<DblMatrix> &matrices, const IntMatrix &mapping, const std::vector<std::string> &labels, const PhyTree &tree) throw (ParameterException);
    //      ~MinSquareTreeCollection();
    void compute(bool KeepTopology, int iter=3, bool quiet=false) throw (RuntimeException);
    bool isInitialized();
    bool isComputed();
    void getTree();
    double getQual() { return this->MST_Qual; };
    PhyTree::TreePtr getPhyTree();
    double getScore();
    double getLogLikelihood();
    std::string newick = "";
    double bestScore;
    double bestLnL;
    PhyTree::TreePtr bestTree;
};

#endif
