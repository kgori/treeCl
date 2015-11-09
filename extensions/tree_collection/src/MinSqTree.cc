#include "MinSqTree.h"
#include "newick.h"

#include <cmath>
#include <cstdio>
#include <cassert>
#include <string>
#include "Eigen/Core"
#include "Eigen/LU"
#include <vector>

#include <cfloat>

using namespace std;

#define dinv(i)                   ((i) == 0 ? 0.0 : 1.0/(i))
#define isNA(k, i)                 (aMap(k,i) == -1)
#define ReallyLessCollection(a, b) (((b)-(a)) > max(1.0e-8*(abs(b)+abs(a)),1.0e-8))
#define ReallyLessGrob(a, b)       (((b)-(a)) > max(1.0e-6*(abs(b)+abs(a)),1.0e-6))
#define nSpecies(k)               (aDistVar[k].rows())
#define aDIST(k, i, j)              (aMap(k,i) < aMap(k,j) ? (aDistVar[k](aMap(k,i)-1,aMap(k,j)-1)) : (aDistVar[k](aMap(k,j)-1,aMap(k,i)-1)))
#define aVAR(k, i, j)               (dinv(aMap(k,i) < aMap(k,j) ? (aDistVar[k](aMap(k,j)-1,aMap(k,i)-1)) : (aDistVar[k](aMap(k,i)-1,aMap(k,j)-1))))
#define VAR(k,i,j)                (aMap(k,i) < aMap(k,j) ? (aDistVar[k](aMap(k,j)-1,aMap(k,i)-1)) : (aDistVar[k](aMap(k,i)-1,aMap(k,j)-1)))

/* if we want to replace the macro by inline function, will need to solve
 * the scope problem of aMap and aDistVar
 inline double aVAR(int k, int i, int j) {
	int i1 = aMap(k,i);
	int j1 = aMap(k,j);
	double var = i1 < j1 ?
 aDistVar[k](j1-1,i1-1)
 : aDistVar[k](i1-1,j1-1)
	return dinv(var);
 }

 inline double aDIST(int k, int i, int j) {
	int i1 = aMap(k,i);
	int j1 = aMap(k,j);
	double dist = i1 < j1 ?
 aDistVar[k](i1-1,j1-1)
 : aDistVar[k](j1-1,i1-1)
	return dist;
 }
 */

/* solve restricted system and return norm of the errors */
void MinSquareTreeCollection::SolveRestricted(const DblMatrix &Z, const DblVector &B, DblVector &L, int n, int restr) {
    DblMatrix ZZ = Z;
    DblVector BB = B;
    int i, j;
    for (i = 0; i < n; i++) {
        if ((restr & (1 << i)) != 0) {
            double maxcol = {0};
            for (j = 0; j < n; j++) {
                if (abs(Z(j, i)) > maxcol) maxcol = abs(Z(j, i));
                ZZ(i, j) = 0;
            }
            ZZ(i, i) = 16 * maxcol;
            BB[i] = MinLen * 16 * maxcol;
        }
        //gausselimi(ZZ,BB,L,n);
        L = ZZ.lu().solve(BB);
    }
}

/* Does brute force search on the possible subset restrictions of restr.
 It is assumed that for restr we have a feasible solution */
int MinSquareTreeCollection::BestRestrictionRemoval(const DblMatrix &Z, const DblVector &B, DblVector &L, int n, int restr) {
    double bestnorm, eps, norm, max;
    int bestrestr, i, k;
    DblVector L2(n);
    DblVector bestL(n);

    bestnorm = DBL_MAX;
    bestrestr = restr;
    for (k = 0; k < n; k++)
        if ((restr & (1 << k)) != 0) {
            SolveRestricted(Z, B, L2, n, restr - (1 << k));
            for (max = i = 0; i < n; i++) if (abs(L2[i]) > max) max = abs(L2[i]);
            eps = 50 * DBL_EPSILON * max;
            for (i = 0; i < n; i++) if (L2[i] < MinLen - eps) break;
            if (i < n) /* removing restriction was not good */ continue;

            bestrestr = BestRestrictionRemoval(Z, B, L2, n, restr - (1 << k));
            if (bestrestr == restr) continue;
            /* The norm of the error is btb - x*Atb (for Ax ~ b) */
            /* we only compute -x*Atb to select the smallest error */
            /* In here Z=AtA, B=Atb, L2=x */
            for (norm = i = 0; i < n; i++) norm -= L2[i] * B[i];
            if (norm < bestnorm) {
                bestnorm = norm;
                bestrestr = restr - (1 << k);
                bestL = L2;
            }
        }
    if (bestrestr != restr) L = bestL;
    return (bestrestr);
}


double MinSquareTreeCollection::FourSubtree(double Wab, double Mab, double Wac, double Mac,
                                            double Wad, double Mad, double Wbc, double Mbc,
                                            double Wbd, double Mbd, double Wcd, double Mcd,
                                            DblVector &L) {
    double eps, max;
    DblMatrix Z1(5, 5);
    DblVector B1(5), L2(5);
    int i, it, restr;
    Z1 <<
    Wab + Wac + Wad,
    Wab,
    Wac + Wad,
    Wac,
    Wad,
    Wab,
    Wab + Wbc + Wbd,
    Wbc + Wbd,
    Wbc,
    Wbd,
    Wac + Wad,
    Wbc + Wbd,
    Wac + Wad + Wbc + Wbd,
    Wac + Wbc,
    Wad + Wbd,
    Wac,
    Wbc,
    Wac + Wbc,
    Wac + Wbc + Wcd,
    Wcd,
    Wad,
    Wbd,
    Wad + Wbd,
    Wcd,
    Wad + Wbd + Wcd;

    B1 <<
    Mab + Mac + Mad,
    Mab + Mbc + Mbd,
    Mac + Mad + Mbc + Mbd,
    Mac + Mbc + Mcd,
    Mad + Mbd + Mcd;

    for (restr = it = 0; it < 5; it++) {
    restart:
        SolveRestricted(Z1, B1, L, 5, restr);
        for (max = i = 0; i < 5; i++) if (abs(L[i]) > max) max = abs(L[i]);
        eps = 50 * DBL_EPSILON * max;
        for (i = 0; i < 5; i++)
            if (L[i] < MinLen - eps) {
                assert((restr & (1 << i)) == 0);
                restr |= 1 << i;
                goto restart;
            }
    }
    /* if there are two or more restrictions to MinLen, try all possible
     combinations of those restrictions */
    if (restr != 0 && (restr & (restr - 1)) != 0) {
        if (restr != BestRestrictionRemoval(Z1, B1, L2, 5, restr)) {
            for (i = 0; i < 5; i++) L[i] = L2[i];
        }
    }
    for (max = i = 0; i < 5; i++) if (abs(L[i]) > max) max = abs(L[i]);
    eps = 50 * DBL_EPSILON * max;
    for (i = 0; i < 5; i++) assert(L[i] >= MinLen - eps);

    {
        double t1, t4, t7, t10, t13, t20;
        t1 = L[2] + L[3] + L[0];
        t4 = L[2] + L[4] + L[0];
        t7 = L[3] + L[2] + L[1];
        t10 = L[4] + L[2] + L[1];
        t13 = L[0] + L[1];
        t20 = L[3] + L[4];
        return (t1 * t1 * Wac + t4 * t4 * Wad + t7 * t7 * Wbc + t10 * t10 * Wbd + t13 * t13 * Wab +
                t20 * t20 * Wcd - 2 * (t1 * Mac + t4 * Mad + t7 * Mbc + t10 * Mbd + t13 * Mab + Mcd * t20));
    }
}


/* The global arrays ConShortestPath[leaf] and ShortestLabel[leaf] are
 set by MS_ShortestPath  */

void MinSquareTreeCollection::printTreeC() {
    for (int k = 0; k < NT; k++) {
        printf("OG(%d)[", k);
        for (int i = 0; i < 2 * ne - 3; i++) {
            if (EdgeC[i].len[k] == DBL_MAX) {
                printf("[e:%d f:%d t:%d l:MAX]", i, EdgeC[i].From, EdgeC[i].To);
            } else {
                printf("[e:%d f:%d t:%d l:%f]", i, EdgeC[i].From, EdgeC[i].To,
                       EdgeC[i].len[k]);
            }
        }
        printf("]\n");
    }
}


void MinSquareTreeCollection::printTree() {
    for (int i = 0; i < 2 * ne - 3; i++)
        printf("[e:%d f:%d t:%d l:%f]", i, Edge[i].From, Edge[i].To,
               Edge[i].len);
    printf("\n");
}


/* this was adapted to store the tmp length in internal nodes */
void MinSquareTreeCollection::MS_ShortestPathCollection(int from, int ExcludedEdge, int label)
/* requires inc[][] vector */
{
    int dest;
    if (from < ne) {
        ShortestLabel[from] = label;
        return;
    }

    for (int j = 0; j < 3; j++) {
        int e1 = inc(from - ne, j);
        if (e1 != ExcludedEdge) {
            if (EdgeC[e1].From == from) {
                dest = EdgeC[e1].To;
            }
            else {
                dest = EdgeC[e1].From;
            }
            for (int k = 0; k < NT; k++) {
                ConShortestPathC(dest, k) = ConShortestPathC(from, k)
                + abs(EdgeC[e1].len[k]);
            }
            MS_ShortestPathCollection(dest, e1, label);
        }
    }
    return;
}

void MinSquareTreeCollection::MS_ShortestPathOne(int from, int ExcludedEdge, int k, int label)
/* requires inc[][] vector */
{
    int orig, dest;
    if (from < ne) {
        ShortestLabel[from] = label;
        return;
    }

    for (int j = 0; j < 3; j++) {
        int e1 = inc(from - ne, j);
        if (e1 != ExcludedEdge) {
            if (EdgeC[e1].From == from) {
                dest = EdgeC[e1].To;
                orig = EdgeC[e1].From;
            }
            else {
                dest = EdgeC[e1].From;
                orig = EdgeC[e1].To;
            }
            ConShortestPathC(dest, k) = ConShortestPathC(orig, k)
            + abs(EdgeC[e1].len[k]);

            MS_ShortestPathOne(dest, e1, k, label);
        }
    }
    return;
}

/* computes the sum of log-likelihoods of the distance errors as given in Edge
*  assumes normally-distributed errors */
double MinSquareTreeCollection::LogLikelihoodFitCollection()
/* requires inc[][] vector */
{
   int e1, e2, e3;
   double error, r, variance, l;
   r = 0;
   for(e3=0; e3 < 2*ne-3; e3++ ) {
      if(    EdgeC[e3].From < ne-1 ){ e1 = EdgeC[e3].From; e2 = EdgeC[e3].To;}
      else if( EdgeC[e3].To < ne-1 ){ e1 = EdgeC[e3].To;   e2 = EdgeC[e3].From;}
      else continue;

      for(int k=0; k<NT; k++) {
         ConShortestPathC(e2,k) = abs(EdgeC[e3].len[k]);
      }
      MS_ShortestPathCollection( e2, e3, 0 );

      for(int k=0; k<NT; k++) {
         if( isNA(k,e1) )
            continue;
         for(int j=e1+1; j<ne; j++ ) {
            if( isNA(k,j) )
               continue;
            error = ConShortestPathC(j,k) - aDIST(k,e1,j);
            variance = VAR(k,e1,j);
            l = -0.5*(log(2*M_PI*variance) + error*error/variance);
            r += l;
         }
      }
   }
   return(r);
}

/* computes the sum of squares of the distance errors as given in Edge */
double MinSquareTreeCollection::DistanceFitCollection()
/* requires inc[][] vector */
{
    int e1, e2, e3;
    double r, t, t2;
    r = 0;
    for (int k = 0; k < NT; k++) {
        tmpA[k] = 0;
    }
    for (e3 = 0; e3 < 2 * ne - 3; e3++) {
        if (EdgeC[e3].From < ne - 1) {
            e1 = EdgeC[e3].From;
            e2 = EdgeC[e3].To;
        }
        else if (EdgeC[e3].To < ne - 1) {
            e1 = EdgeC[e3].To;
            e2 = EdgeC[e3].From;
        }
        else continue;

        for (int k = 0; k < NT; k++) {
            ConShortestPathC(e2, k) = abs(EdgeC[e3].len[k]);
        }
        MS_ShortestPathCollection(e2, e3, 0);

        for (int k = 0; k < NT; k++) {
            if (isNA(k, e1))
                continue;
            for (int j = e1 + 1; j < ne; j++) {
                if (isNA(k, j))
                    continue;
                t = ConShortestPathC(j, k) - aDIST(k, e1, j);
                t2 = t * t * aVAR(k, e1, j);
                r += t2;
                tmpA[k] += t2;
            }
        }
    }
    return (r);
}

void MinSquareTreeCollection::IncidencesC() {
    for (int v1 = ne - 3; v1 >= 0; v1--)
        for (int j = 0; j < 3; j++)
            inc(v1, j) = (-1);
    for (int e0 = 2 * ne - 4; e0 >= 0; e0--) {
        if (EdgeC[e0].From >= ne) {
            int i0 = EdgeC[e0].From - ne;
            if (inc(i0, 0) == (-1)) inc(i0, 0) = e0;
            else if (inc(i0, 1) == (-1)) inc(i0, 1) = e0;
            else inc(i0, 2) = e0;
        }
        if (EdgeC[e0].To >= ne) {
            int i0 = EdgeC[e0].To - ne;
            if (inc(i0, 0) == (-1)) inc(i0, 0) = e0;
            else if (inc(i0, 1) == (-1)) inc(i0, 1) = e0;
            else inc(i0, 2) = e0;
        }
    }
    return;
}

/*

 Code to fit a branches of subtrees of orthologous groups. All branches
 labeled with negative length in the len-field are optimized.

 Main procedure: FitLabeledEdgesC
 Sub routines: getSons, LabelNonExistEdgesR, LabelNonExistEdges,
 FindQuartet, CountOrLabelPath, FitQuartet

 */

void MinSquareTreeCollection::getSons(int e0, int n, int *e1, int *e2) {
    int k;
    k = n - ne;
    if (inc(k, 0) == e0) {
        *e1 = inc(k, 1);
        *e2 = inc(k, 2);
    }
    else if (inc(k, 1) == e0) {
        *e1 = inc(k, 0);
        *e2 = inc(k, 2);
    }
    else {
        *e1 = inc(k, 0);
        *e2 = inc(k, 1);
    }
    return;
}

void MinSquareTreeCollection::LabelNonExistEdgesR(int e0, int n_papa) {
    int n, e1, e2;

    n = (EdgeC[e0].From != n_papa) ? EdgeC[e0].From : EdgeC[e0].To;

    if (n < ne) {
        for (int k = 0; k < NT; k++) {
            if (isNA(k, n)) EdgeC[e0].len[k] = DBL_MAX;
        }
        return;
    }
    /* remove DBL_MAX flag from previous run */
    for (int k = 0; k < NT; k++) {
        if (EdgeC[e0].len[k] == DBL_MAX)
            EdgeC[e0].len[k] = -1;
    }

    getSons(e0, n, &e1, &e2);

    LabelNonExistEdgesR(e1, n);
    LabelNonExistEdgesR(e2, n);

    for (int k = 0; k < NT; k++) {
        if ((EdgeC[e1].len[k] == DBL_MAX) && (EdgeC[e2].len[k] == DBL_MAX))
            EdgeC[e0].len[k] = DBL_MAX;
    }
    return;
}

void MinSquareTreeCollection::Fix001Case(int er, int e0, int n_papa) {
    int e1, e2, n;

    n = (EdgeC[e0].From != n_papa) ? EdgeC[e0].From : EdgeC[e0].To;
    if (n < ne) return;

    getSons(e0, n, &e1, &e2);

    for (int k = 0; k < NT; k++) {
        if (EdgeC[er].len[k] == -DBL_MAX) {
            if (e0 != er) EdgeC[e0].len[k] = DBL_MAX;
            if (EdgeC[e1].len[k] != DBL_MAX && EdgeC[e2].len[k] != DBL_MAX)
                EdgeC[er].len[k] = DBL_MAX;
        }
    }

    Fix001Case(er, e1, n);
    Fix001Case(er, e2, n);

    return;
}

void MinSquareTreeCollection::LabelNonExistEdges() {
    int e0, e1, e2, n = 0, n2 = 0;

    for (e0 = 0; e0 <= 2 * ne - 4; e0++) {
        n = EdgeC[e0].From;
        n2 = EdgeC[e0].To;
        if (n < ne) break;
        n = EdgeC[e0].To;
        n2 = EdgeC[e0].From;
        if (n < ne) break;
    }


    LabelNonExistEdgesR(e0, n);

    getSons(e0, n2, &e1, &e2);
    for (int k = 0; k < NT; k++) {
        if (isNA(k, n)) {
            EdgeC[e0].len[k] = DBL_MAX; /* because not labeled by
                                         LabelNonExistEdgesR */
            if ((EdgeC[e1].len[k] == DBL_MAX)
                || (EdgeC[e2].len[k] == DBL_MAX)) {
                EdgeC[e0].len[k] = -DBL_MAX;
            }
        }
    }

    Fix001Case(e0, e0, n);
    for (int k = 0; k < NT; k++)
        if (EdgeC[e0].len[k] == -DBL_MAX)
            throw RuntimeException("LabelNonExistEdges -- inconsistency");

    return;
}

void MinSquareTreeCollection::FindQuartet(int i, int e0, int n_papa, int *at, int *et, int *x, int *ia,
                                          int *ix, int *a2, int *ee2, int *ia2, int *ex2) {
    int n_next, e1, e2, n1, n2;

    n_next = (EdgeC[e0].To == n_papa) ? EdgeC[e0].From : EdgeC[e0].To;

    /* Leaf */
    if (n_next < ne) {
        (*ia)++;
        at[(*ia) - 1] = n_next;
        et[(*ia) - 1] = e0;
        /* the following is to track which at[.] are cherries */
        if (*ix == 2 && (*ia2) < 2) {
            (*ia2)++;
            a2[(*ia2) - 1] = n_next;
            ee2[(*ia2) - 1] = e0;
        }
        return;
    }

    getSons(e0, n_next, &e1, &e2);
    n1 = (EdgeC[e1].To == n_next) ? EdgeC[e1].From : EdgeC[e1].To;
    n2 = (EdgeC[e2].To == n_next) ? EdgeC[e2].From : EdgeC[e2].To;

    /* Path */
    if (EdgeC[e1].len[i] == DBL_MAX) {
        FindQuartet(i, e2, n_next, at, et, x, ia, ix, a2, ee2, ia2, ex2);
        return;
    } else if (EdgeC[e2].len[i] == DBL_MAX) {
        FindQuartet(i, e1, n_next, at, et, x, ia, ix, a2, ee2, ia2, ex2);
        return;
    }

    /* Split */
    if (*ix < 2) {
        (*ix)++;
        x[(*ix) - 1] = n_next;
        if (*ix == 2) *ex2 = e0;
        FindQuartet(i, e1, n_next, at, et, x, ia, ix, a2, ee2, ia2, ex2);
        FindQuartet(i, e2, n_next, at, et, x, ia, ix, a2, ee2, ia2, ex2);

        return;
    } else {
        (*ia)++;
        at[(*ia) - 1] = n_next;
        et[(*ia) - 1] = e0;
        if (*ix == 2 && (*ia2) < 2) {
            (*ia2)++;
            a2[(*ia2) - 1] = n_next;
            ee2[(*ia2) - 1] = e0;
        }
        return;
    }
}

void MinSquareTreeCollection::FindSplit(int i, int e0, int n_papa, int *x, int *ix, int *at, int *ia,
                                        int *et, int *ex) {
    int n_next, e1, e2, n1, n2;

    n_next = (EdgeC[e0].To == n_papa) ? EdgeC[e0].From : EdgeC[e0].To;
    if (n_next < ne) {
        (*ia)++;
        at[(*ia) - 1] = n_next;
        et[(*ia) - 1] = e0;
        return;
    }

    getSons(e0, n_next, &e1, &e2);
    n1 = (EdgeC[e1].To == n_next) ? EdgeC[e1].From : EdgeC[e1].To;
    n2 = (EdgeC[e2].To == n_next) ? EdgeC[e2].From : EdgeC[e2].To;

    if (EdgeC[e1].len[i] == DBL_MAX) {
        FindSplit(i, e2, n_next, x, ix, at, ia, et, ex);
        return;
    } else if (EdgeC[e2].len[i] == DBL_MAX) {
        FindSplit(i, e1, n_next, x, ix, at, ia, et, ex);
        return;
    }

    /* Split */
    if (*ix < 1) { /* inner split */
        (*ix)++;
        x[(*ix) - 1] = n_next;
        *ex = e0;
        FindSplit(i, e1, n_next, x, ix, at, ia, et, ex);
        FindSplit(i, e2, n_next, x, ix, at, ia, et, ex);
        return;
    } else { /* outer split */
        (*ia)++;
        at[(*ia) - 1] = n_next;
        et[(*ia) - 1] = e0;
        return;
    }
}

int MinSquareTreeCollection::CountOrLabelPath(int i, int x2, int x1, int ex2, double dL, int *ise0) {
    int e1, e2, np = 0, ep, l;
    int isToBeFitted = 0;
    l = 0;

    if (EdgeC[ex2].From == x2) {np = EdgeC[ex2].To;}
    else if (EdgeC[ex2].To == x2) {np = EdgeC[ex2].From;}
    else {throw RuntimeException("CountOrLabelPath -- inconsistent input");}
    ep = ex2;
    if (ep == (*ise0)) (*ise0) = -1;
    if (dL != -1) EdgeC[ep].len[i] = dL;
    l++;

    while (np != x1) {
        if (np < ne) {
            /*            printf("x2 %d x1 %d ex2 %d\n", x2, x1, ex2); */
            throw RuntimeException("CountOrLabelPath -- bug: reached a leaf");
        }
        getSons(ep, np, &e1, &e2);
        if (EdgeC[e1].len[i] == DBL_MAX) {ep = e2;} else {ep = e1;}
        if (ep == (*ise0)) (*ise0) = -1;
        if (dL != -1) EdgeC[ep].len[i] = dL;
        l++;
        /* negative l means that there is something to be fitted */
        if (EdgeC[e1].len[i] < 0) {
            isToBeFitted = 1;
        }
        if (EdgeC[ep].From == np) {np = EdgeC[ep].To;}
        else {np = EdgeC[ep].From;}
    }

    return (isToBeFitted ? -l : l);
}

void MinSquareTreeCollection::ThreeOptimSubset(double wab, double mab, double wac, double mac, double wbc,
                                               double mbc, double *T) {
    double min, tmp, t1, t2, t3, ab, ac, bc;
    ab = mab / wab;
    ac = mac / wac;
    bc = mbc / wbc;
    T[0] = (ab + ac - bc) / 2;
    T[1] = (ab - ac + bc) / 2;
    T[2] = (-ab + ac + bc) / 2;
    if (T[0] < 0 || T[1] < 0 || T[2] < 0) {
        // 0 0 0
        T[0] = T[1] = T[2] = 0;
        min = pow(ab, 2) + pow(ac, 2) + pow(bc, 2);
        // 0 0 1
        t1 = (mac + mbc) / (wac + wbc);
        tmp = pow(ab, 2) + pow(ac - t1, 2) + pow(bc - t1, 2);
        if (tmp < min && t1 > 0) {
            T[2] = t1;
            min = tmp;
        }
        // 0 1 0
        t1 = (mab + mbc) / (wab + wbc);
        tmp = pow(ab - t1, 2) + pow(ac, 2) + pow(bc - t1, 2);
        if (tmp < min && t1 > 0) {
            T[0] = T[2] = 0;
            T[1] = t1;
            min = tmp;
        }
        // 1 0 0
        t1 = (mab + mac) / (wab + wac);
        tmp = pow(ab - t1, 2) + pow(ac - t1, 2) + pow(bc, 2);
        if (tmp < min && t1 > 0) {
            T[1] = T[2] = 0;
            T[0] = t1;
            min = tmp;
        }

        t3 = pow(wab * wbc, 2) + pow(wac * wab, 2) + pow(wac * wbc, 2);
        // 0 1 1
        t1 = (-wbc * wbc * wac * wac * ac + wab * wab * ab * wac * wac + wac
              * wac * wbc * wbc * bc + wbc * wbc * wab * wab * ab) / t3;
        t2 = -(-wab * wab * wac * wac * ac - wbc * wbc * wac * wac * ac + wbc
               * wbc * wab * wab * ab - wab * wab * wbc * wbc * bc) / t3;


        tmp = pow(ab - t1, 2) + pow(ac - t2, 2) + pow(bc - t1 - t2, 2);
        if (tmp < min && t1 > 0 && t2 > 0) {
            T[0] = 0;
            T[1] = t1;
            T[2] = t2;
            min = tmp;
        }
        // 1 0 1
        t1 = (wab * wab * ab * wac * wac - wac * wac * wbc * wbc * bc + wbc *
              wbc * wab * wab * ab + wbc * wbc * wac * wac * ac) / t3;
        t2 = -(wab * wab * ab * wac * wac - wab * wab * wac * wac * ac - wab
               * wab * wbc * wbc * bc - wac * wac * wbc * wbc * bc) / t3;
        tmp = pow(ab - t1, 2) + pow(ac - t1 - t2, 2) + pow(bc - t2, 2);
        if (tmp < min && t1 > 0 && t2 > 0) {
            T[1] = 0;
            T[0] = t1;
            T[2] = t2;
            min = tmp;
        }
        // 1 1 0
        t1 = (-wab * wab * wbc * wbc * bc + wab * wab * wac * wac * ac + wbc * wbc
              * wab * wab * ab + wbc * wbc * wac * wac * ac) / t3;
        t2 = (wab * wab * wbc * wbc * bc - wab * wab * wac * wac * ac + wac *
              wac * wbc * wbc * bc + wab * wab * ab * wac * wac) / t3;
        tmp = pow(ab - t1 - t2, 2) + pow(ac - t1, 2) + pow(bc - t2, 2);
        if (tmp < min && t1 > 0 && t2 > 0) {
            T[2] = 0;
            T[0] = t1;
            T[1] = t2;
            min = tmp;
        }
    }
}

int MinSquareTreeCollection::CountOrLabelPathTriplet(int i, int A, int x0, int eA, double dL) {
    int e1, e2, np = 0, ep, l;
    l = 1;

    if (EdgeC[eA].From == A) {np = EdgeC[eA].To;}
    else if (EdgeC[eA].To == A) {np = EdgeC[eA].From;}
    else {throw RuntimeException("CountOrLabelPathTriplet -- inconsistent input");}
    ep = eA;

    EdgeC[ep].len[i] = dL;
    while (np != x0) {
        getSons(ep, np, &e1, &e2);
        if (EdgeC[e1].len[i] == DBL_MAX) {
            ep = e2;
        } else if (EdgeC[e2].len[i] == DBL_MAX) {
            ep = e1;
        } else {
            break;
        }
        EdgeC[ep].len[i] = dL;
        l++;
        if (EdgeC[ep].From == np) {np = EdgeC[ep].To;}
        else {np = EdgeC[ep].From;}
    }

    return (l);
}

#define a 0
#define b 1
#define c 2
#define d 4
#define e 7
#define Mab M[a+b]
#define Mac M[a+c]
#define Mad M[a+d]
#define Mae M[a+e]
#define Mbc M[b+c]
#define Mbd M[b+d]
#define Mbe M[b+e]
#define Mcd M[c+d]
#define Mce M[c+e]
#define Mde M[d+e]
#define Wab W[a+b]
#define Wac W[a+c]
#define Wad W[a+d]
#define Wae W[a+e]
#define Wbc W[b+c]
#define Wbd W[b+d]
#define Wbe W[b+e]
#define Wcd W[c+d]
#define Wce W[c+e]
#define Wde W[d+e]


void MinSquareTreeCollection::FitTriplet(int k, int A, int B, int C, int eA, int eB, int eC, int x0) {
    double W[7], M[7], T[3], dd, Wei;

    ConShortestPathC(A, k) = 0;
    ConShortestPathC(B, k) = 0;
    ConShortestPathC(C, k) = 0;
    MS_ShortestPathOne(A, eA, k, a);
    MS_ShortestPathOne(B, eB, k, b);
    MS_ShortestPathOne(C, eC, k, c);

    for (int j = 0; j <= 6; j++) {W[j] = M[j] = 0.0;}

    for (int si = 0; si < ne; si++) {
        if (!isNA(k, si)) {
            for (int sj = si + 1; sj < ne; sj++) {
                if (ShortestLabel[si] != ShortestLabel[sj]) {
                    if (!isNA(k, sj)) {
                        int i = ShortestLabel[si] + ShortestLabel[sj];
                        double v = aVAR(k, si, sj);
                        if (v == 0) {
                            Wei = 1.0e-10;
                            dd = 300.0;
                        }
                        else {
                            Wei = v;
                            dd = aDIST(k, si, sj);
                        }
                        W[i] += Wei;
                        M[i] += Wei * (dd - ConShortestPathC(si, k)
                                       - ConShortestPathC(sj, k));
                    }
                }
            }
        }
    }

    ThreeOptimSubset(W[a + b], M[a + b], W[a + c], M[a + c], W[b + c], M[b + c], T);
    int l = CountOrLabelPathTriplet(k, A, x0, eA, 0);
    CountOrLabelPathTriplet(k, A, x0, eA, T[0] / l);

    return;
}


void MinSquareTreeCollection::FitQuartet(int k, int A, int B, int C, int D, int eA, int eB, int eC,
                                         int eD, int x1, int x2, int ex2, int qOk, int e0, int allEdges) {
    int l, e0f;
    double W[7], M[7], dd, Wei;
    DblVector L(5);

    ConShortestPathC(A, k) = 0;
    ConShortestPathC(B, k) = 0;
    ConShortestPathC(C, k) = 0;
    ConShortestPathC(D, k) = 0;
    MS_ShortestPathOne(A, eA, k, a);
    MS_ShortestPathOne(B, eB, k, b);
    MS_ShortestPathOne(C, eC, k, c);
    MS_ShortestPathOne(D, eD, k, d);

    for (int j = 0; j <= 6; j++) {W[j] = M[j] = 0.0;}

    for (int si = 0; si < ne; si++) {
        if (!isNA(k, si)) {
            for (int sj = si + 1; sj < ne; sj++) {
                if (ShortestLabel[si] != ShortestLabel[sj]) {
                    if (!isNA(k, sj)) {
                        int i = ShortestLabel[si] + ShortestLabel[sj];
                        double v = aVAR(k, si, sj);
                        if (v == 0) {
                            Wei = 1.0e-10;
                            dd = 300.0;
                        }
                        else {
                            Wei = v;
                            dd = aDIST(k, si, sj);
                        }
                        W[i] += Wei;
                        M[i] += Wei * (dd - ConShortestPathC(si, k)
                                       - ConShortestPathC(sj, k));
                    }
                }
            }
        }
    }

    FourSubtree(W[a + b], M[a + b],
                W[a + c], M[a + c],
                W[a + d], M[a + d],
                W[b + c], M[b + c],
                W[b + d], M[b + d],
                W[c + d], M[c + d], L);


    /* middle branch */
    l = abs(CountOrLabelPath(k, x2, x1, ex2, 0, &e0f));
    CountOrLabelPath(k, x2, x1, ex2, L[2] / l, &e0f);
    /**** WARNING: the following will never be executed !!!!! ******/
    if (!allEdges) {
        /*           throw RuntimeException("FitQuartet -- should never be here"); */
        e0f = e0;
        l = CountOrLabelPath(k, A, x1, eA, -1, &e0f);
        if (l > 0 || e0f == -1)
            CountOrLabelPath(k, A, x1, eA, L[0] / abs(l), &e0f);
        else
            CountOrLabelPath(k, A, x1, eA, L[0] / l, &e0f);

        e0f = e0;
        l = CountOrLabelPath(k, B, x1, eB, -1, &e0f);
        if (l > 0 || e0f == -1)
            CountOrLabelPath(k, B, x1, eB, L[1] / abs(l), &e0f);
        else
            CountOrLabelPath(k, B, x1, eB, L[1] / l, &e0f);

        e0f = e0;
        l = CountOrLabelPath(k, C, x2, eC, -1, &e0f);
        if (l > 0 || e0f == -1)
            CountOrLabelPath(k, C, x2, eC, L[3] / abs(l), &e0f);
        else
            CountOrLabelPath(k, C, x2, eC, L[3] / l, &e0f);

        e0f = e0;
        l = CountOrLabelPath(k, D, x2, eD, -1, &e0f);
        if (l > 0 || e0f == -1)
            CountOrLabelPath(k, D, x2, eD, L[4] / abs(l), &e0f);
        else
            CountOrLabelPath(k, D, x2, eD, L[4] / l, &e0f);
    }
    if (allEdges) {
        l = abs(CountOrLabelPath(k, A, x1, eA, -1, &e0f));
        CountOrLabelPath(k, A, x1, eA, L[0] / l, &e0f);

        l = abs(CountOrLabelPath(k, B, x1, eB, -1, &e0f));
        CountOrLabelPath(k, B, x1, eB, L[1] / l, &e0f);

        l = abs(CountOrLabelPath(k, C, x2, eC, -1, &e0f));
        CountOrLabelPath(k, C, x2, eC, L[3] / l, &e0f);

        l = abs(CountOrLabelPath(k, D, x2, eD, -1, &e0f));
        CountOrLabelPath(k, D, x2, eD, L[4] / l, &e0f);
    }

    return;
}

/*



 A = at[i]                  C = at[j]
 \ et[i]                 / et[j]
 o                     o
 '                     '
 '                     '
 o                     o
 \     3             /
 x[1] o----o.....o------o x[2]
 /               ex2 \
 o                     o
 '                     '
 / et[l]                 \ et[k]
 B = at[l]                  D = at[k]

 */

void MinSquareTreeCollection::GotoLeaf(int i, int e0, int n_papa, int *eL, int *L) {
    int np = 0, ep, e1, e2;

    if (n_papa < ne) {
        (*eL) = e0;
        (*L) = n_papa;
        return;
    }

    if (EdgeC[e0].From == n_papa) {np = EdgeC[e0].To;}
    else if (EdgeC[e0].To == n_papa) {np = EdgeC[e0].From;}
    else {throw RuntimeException("GotoLeaf -- inconsistent input");}

    if (np < ne) {
        (*eL) = e0;
        (*L) = np;
        return;
    }

    getSons(e0, np, &e1, &e2);
    if (EdgeC[e1].len[i] == DBL_MAX) {ep = e2;} else {ep = e1;}
    GotoLeaf(i, ep, np, eL, L);
}

#if 0
void MinSquareTreeCollection::FitEdge_SAVE(int i, int e0, int allEdges)
{
    int j, at[4], et[4], x[2], ia, ix, a2[2], ee2[2], ex2, ia2,
    A, B, C, D, eA, eB, eC, eD, qOk, ex;
    int k, trip;

    qOk = 0;
    ia = ix = 0;
    trip = 0;
    FindSplit(i, e0, EdgeC[e0].To, x, &ix, at, &ia, et, &ex);
    if(ix == 1) {
        C  = at[0];  D = at[1];  x[1] = x[0];
        eC = et[0]; eD = et[1];  ex2 = ex;
        ia = ix = 0;
        FindSplit(i, e0, EdgeC[e0].From, x, &ix, at, &ia, et, &ex);
        if(ix == 1) {
            qOk = 1;
            A  = at[0]; B  = at[1]; /* x[0] is set in FindSplit */
            eA = et[0]; eB = et[1];
        }
    } else {
        trip = 1;
        ia = ix = 0;
        FindSplit(i, e0, EdgeC[e0].From, x, &ix, at, &ia, et, &ex);
        if(ix == 1) {
            A  = at[0]; B  = at[1]; /* x[0] is set in FindSplit */
            eA = et[0]; eB = et[1];
        } else {
            throw RuntimeException("FitEdge -- internal inconsistency");
        }
    }

    if(qOk == 0) {
        if(trip == 0) {
            GotoLeaf(i, e0, EdgeC[e0].From, &eA, &A);
            FitTriplet(i, A, C, D, eA, eC, eD, x[1]);
            /* old: FitTriplet(e0, EdgeC[e0].To,   C, D, eC, eD, x[1]); */
        } else {
            GotoLeaf(i, e0, EdgeC[e0].To, &eC, &C);
            FitTriplet(i, C, A, B, eC, eA, eB, x[0]);
            /* old: FitTriplet(e0, EdgeC[e0].From, A, B, eA, eB, x[0]); */
        }
    } else {
        FitQuartet(i, A, B, C, D, eA, eB, eC, eD, x[1-1], x[2-1], ex2, qOk, e0,
                   allEdges);
    }

    return;
}
#endif

void MinSquareTreeCollection::FitEdge(int i, int e0, int allEdges) {
    int at[4], et[4], x[2], ia, ix, a2[2], ee2[2], ex2, ia2,
    A, B, C, D, eA = 0, eB = 0, eC, eD, qOk, ex;

    qOk = 0;
    ia = ix = 0;
    FindSplit(i, e0, EdgeC[e0].To, x, &ix, at, &ia, et, &ex);
    if (ix == 1) {
        C = at[0];
        D = at[1];
        x[1] = x[0];
        eC = et[0];
        eD = et[1];
        ex2 = ex;
        ia = ix = 0;
        FindSplit(i, e0, EdgeC[e0].From, x, &ix, at, &ia, et, &ex);
        if (ix == 1) {
            qOk = 1;
            A = at[0];
            B = at[1]; /* x[0] is set in FindSplit */
            eA = et[0];
            eB = et[1];
        }
    }

    if (qOk == 0) {
        ia = ix = ia2 = 0;
        FindQuartet(i, e0, EdgeC[e0].To, at, et, x, &ia, &ix, a2, ee2, &ia2,
                    &ex2);
        if (ia < 4)
            FindQuartet(i, e0, EdgeC[e0].From, at, et, x, &ia, &ix, a2, ee2,
                        &ia2, &ex2);

        /* group the cherries ((A,B),(C,D)) in the at[i] */
        A = B = -1;
        C = a2[1 - 1];
        D = a2[2 - 1];
        eC = ee2[1 - 1];
        eD = ee2[2 - 1];
        for (int j = 0; j < 4; j++) {
            if (at[j] != a2[1 - 1] && at[j] != a2[2 - 1] && A == -1) {
                A = at[j];
                eA = et[j];
            } else if (at[j] != a2[1 - 1] && at[j] != a2[2 - 1] && A != -1) {
                B = at[j];
                eB = et[j];
            }
        }
    }

    FitQuartet(i, A, B, C, D, eA, eB, eC, eD, x[1 - 1], x[2 - 1], ex2, qOk, e0,
               allEdges);

    return;
}

void MinSquareTreeCollection::FitLabeledEdgesC(int allEdges) {
    if (allEdges == 1) {
        IncidencesC();
        LabelNonExistEdges();
        for (int k = 0; k < NT; k++) {
            for (int e0 = 0; e0 <= 2 * ne - 4; e0++) {
                if (EdgeC[e0].len[k] != DBL_MAX) {
                    FitEdge(k, e0, allEdges);
                }
            }
            /*

             for(e0=0; e0<=2*ne-4; e0++) {
             if(EdgeC[e0].len[i] != DBL_MAX) {
             if(EdgeC[e0].From >= ne && EdgeC[e0].To >= ne) {
             FitEdge(i, e0, allEdges);
             }
             }
             }

             for(e0=0; e0<=2*ne-4; e0++) {
             if(EdgeC[e0].len[i] != DBL_MAX) {
             if(EdgeC[e0].From < ne && EdgeC[e0].To < ne) {
             FitEdge(i, e0, allEdges);
             }
             }
             }
             */
        }
    } else {
        LabelNonExistEdges();
        for (int k = 0; k < NT; k++) {
            for (int e0 = 0; e0 <= 2 * ne - 4; e0++) {
                if (EdgeC[e0].len[k] != DBL_MAX && EdgeC[e0].len[k] < 0) {
                    FitEdge(k, e0, allEdges);
                }
            }
        }
    }
    return;
}

int MinSquareTreeCollection::FourOptimCollection(int DoNotModify, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W) {
    int e1, e2, e3, e4, e5, ix, iy, va, vb, vc,
    vd, flag[5], tot, nswap;
    double dd, t1, t2, t3, Wei, tmp; //,lik1,lik2,lik3;
    DblVector Lout = DblVector::Zero(7);

    /* Four subtree case:
     #
     #   a                  c
     #    \                /
     #    1\              /4
     #      \     3      /
     #       o----------o
     #      /x          y\
     #    2/              \5
     #    /                \
     #   b                  d
     #
     #
     # L: [k][5]  M,W: [k][7]
     #                           */

    nswap = 0;
    IncidencesC();
    for (e3 = 2 * ne - 4; e3 >= 0; e3--)
        if (EdgeC[e3].From >= ne && EdgeC[e3].To >= ne) {
            ix = EdgeC[e3].From - ne;
            iy = EdgeC[e3].To - ne;

            if (inc(ix, 0) == e3) {
                e1 = inc(ix, 1);
                e2 = inc(ix, 2);
            }
            else if (inc(ix, 1) == e3) {
                e1 = inc(ix, 0);
                e2 = inc(ix, 2);
            }
            else {
                e1 = inc(ix, 0);
                e2 = inc(ix, 1);
            }

            if (inc(iy, 0) == e3) {
                e4 = inc(iy, 1);
                e5 = inc(iy, 2);
            }
            else if (inc(iy, 1) == e3) {
                e4 = inc(iy, 0);
                e5 = inc(iy, 2);
            }
            else {
                e4 = inc(iy, 0);
                e5 = inc(iy, 1);
            }

            va = EdgeC[e1].From == ix + ne ? EdgeC[e1].To : EdgeC[e1].From;
            vb = EdgeC[e2].From == ix + ne ? EdgeC[e2].To : EdgeC[e2].From;
            vc = EdgeC[e4].From == iy + ne ? EdgeC[e4].To : EdgeC[e4].From;
            vd = EdgeC[e5].From == iy + ne ? EdgeC[e5].To : EdgeC[e5].From;

            for (int k = 0; k < NT; k++) {
                ConShortestPathC(va, k) = 0;
                ConShortestPathC(vb, k) = 0;
                ConShortestPathC(vc, k) = 0;
                ConShortestPathC(vd, k) = 0;
            }
            MS_ShortestPathCollection(va, e1, a);
            MS_ShortestPathCollection(vb, e2, b);
            MS_ShortestPathCollection(vc, e4, c);
            MS_ShortestPathCollection(vd, e5, d);

            for (int k = 0; k < NT; k++) {
                for (int i = 0; i <= 6; i++) {
                    W(k, i) = M(k, i) = 0.0;
                }
                flag[a] = flag[b] = flag[c] = flag[d] = 0;

                for (int si = 0; si < ne; si++) {
                    if (isNA(k, si))
                        continue;
                    flag[ShortestLabel[si]] += 1;
                    for (int sj = si + 1; sj < ne; sj++) {
                        if (ShortestLabel[si] != ShortestLabel[sj]) {
                            if (isNA(k, sj))
                                continue;
                            int i = ShortestLabel[si] + ShortestLabel[sj];
                            double v = aVAR(k, si, sj);
                            if (v == 0) {
                                Wei = 1.0e-10;
                                dd = 300.0;
                            }
                            else {
                                Wei = v;
                                dd = aDIST(k, si, sj);
                            }
                            W(k, i) += Wei;
                            M(k, i) += Wei * (dd - ConShortestPathC(si, k) - ConShortestPathC(sj, k));
                        }
                    }
                }
                /* if we are missing too many pairwise distances, the family does
                 not contribue information    */
                if (flag[a] && flag[b] && flag[c] && flag[d])
                    Skip[k] = 0;
                else if (flag[d] && flag[a] && flag[b])
                    Skip[k] = e5 + 1;
                else if ((flag[c] && flag[a] && flag[b])
                         || (flag[a] && flag[b] && !flag[d])) /*c arbitrary*/
                    Skip[k] = e4 + 1;
                else if (flag[b] && ((flag[c] && flag[d])
                                     || (flag[c] && !flag[a])
                                     || (flag[d] && !flag[a])))
                    Skip[k] = e2 + 1;
                else Skip[k] = e1 + 1;
            }
            tot = 0;
            for (int k = 0; k < NT; k++)
                if (Skip[k] == 0)
                    tot++;

            /* t1 */
            t1 = 0;
            if (tot > 0) {
                EdgeC[e1].len0 = 0;
                EdgeC[e2].len0 = 0;
                EdgeC[e3].len0 = 0;
                EdgeC[e4].len0 = 0;
                EdgeC[e5].len0 = 0;
            }
            for (int k = 0; k < NT; k++) {
                if (Skip[k] > 0)
                    continue;
                t1 += FourSubtree(W(k, a + b), M(k, a + b),
                                  W(k, a + c), M(k, a + c),
                                  W(k, a + d), M(k, a + d),
                                  W(k, b + c), M(k, b + c),
                                  W(k, b + d), M(k, b + d),
                                  W(k, c + d), M(k, c + d),
                                  Lout);
                L.row(k) = Lout;
                EdgeC[e1].len[k] = L(k, 0);
                EdgeC[e2].len[k] = L(k, 1);
                EdgeC[e3].len[k] = L(k, 2);
                EdgeC[e4].len[k] = L(k, 3);
                EdgeC[e5].len[k] = L(k, 4);
                EdgeC[e1].len0 += L(k, 0) / tot;
                EdgeC[e2].len0 += L(k, 1) / tot;
                EdgeC[e3].len0 += L(k, 2) / tot;
                EdgeC[e4].len0 += L(k, 3) / tot;
                EdgeC[e5].len0 += L(k, 4) / tot;

            }

            /* t2: (a,c)(b,d)  */
            if (!DoNotModify) {
                t2 = 0;
                for (int k = 0; k < NT; k++) {
                    if (Skip[k] > 0)
                        continue;
                    t2 += FourSubtree(W(k, a + c), M(k, a + c),
                                      W(k, a + b), M(k, a + b),
                                      W(k, a + d), M(k, a + d),
                                      W(k, b + c), M(k, b + c),
                                      W(k, c + d), M(k, c + d),
                                      W(k, b + d), M(k, b + d),
                                      Lout);
                    L.row(k) = Lout;
                }

                if (ReallyLessCollection(t2, t1)) {
                    t1 = t2;
                    nswap++;
                    for (int k = 0; k < NT; k++) {
                        if (Skip[k] > 0) {
                            /* the following could be a bug: what if skip=e5+1? */
                            EdgeC[Skip[k] - 1].len[k] +=
                            EdgeC[e3].len[k];
                            EdgeC[e3].len[k] = 0;
                            tmp = EdgeC[e2].len[k];
                            EdgeC[e2].len[k] = EdgeC[e4].len[k];
                            EdgeC[e4].len[k] = tmp;
                        }
                        else {
                            EdgeC[e1].len[k] = L(k, 0);
                            EdgeC[e2].len[k] = L(k, 1);
                            EdgeC[e3].len[k] = L(k, 2);
                            EdgeC[e4].len[k] = L(k, 3);
                            EdgeC[e5].len[k] = L(k, 4);
                            EdgeC[e1].len0 += L(k, 0) / tot;
                            EdgeC[e2].len0 += L(k, 1) / tot;
                            EdgeC[e3].len0 += L(k, 2) / tot;
                            EdgeC[e4].len0 += L(k, 3) / tot;
                            EdgeC[e5].len0 += L(k, 4) / tot;
                        }
                    }
                    EdgeC[e1].From = ix + ne;
                    EdgeC[e1].To = va;
                    EdgeC[e2].From = ix + ne;
                    EdgeC[e2].To = vc;
                    EdgeC[e4].From = iy + ne;
                    EdgeC[e4].To = vb;
                    EdgeC[e5].From = iy + ne;
                    EdgeC[e5].To = vd;
                    IncidencesC();
#ifdef DEBUG
                    tmp = DistanceFitCollection();
                    printf("FourOptim: 1, edge %d, %.2f\n",e3,tmp);
                    for(int k=0; k<NT; k++) {
                        ConShortestPathC(va,k) = 0;
                        ConShortestPathC(vb,k) = 0;
                        ConShortestPathC(vc,k) = 0;
                        ConShortestPathC(vd,k) = 0;
                    }
                    MS_ShortestPathCollection( va, e1, a );
                    MS_ShortestPathCollection( vc, e2, b );
                    MS_ShortestPathCollection( vb, e4, c );
                    MS_ShortestPathCollection( vd, e5, d );
                    if (ReallyLessCollection(globmin,tmp)) {
                        throw RuntimeException("min increased");
                    }
                    else {
                        globmin = tmp;
                        for(int k=0; k<NT; k++) {
                            globminA[k] = tmpA[k];
                        }
                    }
#endif

                    /* this is because we need to recompute the indidences */
                    e3++;
                    continue;
                }

                /* t3: (a,d)(b,c) */
                t3 = 0;
                for (int k = 0; k < NT; k++) {
                    if (Skip[k] > 0)
                        continue;
                    t3 += FourSubtree(W(k, a + d), M(k, a + d),
                                      W(k, a + c), M(k, a + c),
                                      W(k, a + b), M(k, a + b),
                                      W(k, c + d), M(k, c + d),
                                      W(k, b + d), M(k, b + d),
                                      W(k, b + c), M(k, b + c),
                                      Lout);
                    L.row(k) = Lout;
                }

                if (ReallyLessCollection(t3, t1)) {
                    t1 = t3;
                    nswap++;
                    for (int k = 0; k < NT; k++) {
                        if (Skip[k] > 0) {
                            EdgeC[Skip[k] - 1].len[k] +=
                            EdgeC[e3].len[k];
                            EdgeC[e3].len[k] = 0;
                            tmp = EdgeC[e2].len[k];
                            EdgeC[e2].len[k] = EdgeC[e5].len[k];
                            EdgeC[e5].len[k] = tmp;
                        }
                        else {
                            EdgeC[e1].len[k] = L(k, 0);
                            EdgeC[e2].len[k] = L(k, 1);
                            EdgeC[e3].len[k] = L(k, 2);
                            EdgeC[e4].len[k] = L(k, 3);
                            EdgeC[e5].len[k] = L(k, 4);
                            EdgeC[e1].len0 += L(k, 0) / tot;
                            EdgeC[e2].len0 += L(k, 1) / tot;
                            EdgeC[e3].len0 += L(k, 2) / tot;
                            EdgeC[e4].len0 += L(k, 3) / tot;
                            EdgeC[e5].len0 += L(k, 4) / tot;
                        }
                    }
                    EdgeC[e1].From = ix + ne;
                    EdgeC[e1].To = va;
                    EdgeC[e2].From = ix + ne;
                    EdgeC[e2].To = vd;
                    EdgeC[e4].From = iy + ne;
                    EdgeC[e4].To = vc;
                    EdgeC[e5].From = iy + ne;
                    EdgeC[e5].To = vb;
                    IncidencesC();
#ifdef DEBUG
                    tmp = DistanceFitCollection();
                    printf("FourOptim: 2, edge %d %.3f\n",e3,tmp);

                    for(int k=0; k<NT; k++) {
                        ConShortestPathC(va,k) = 0;
                        ConShortestPathC(vb,k) = 0;
                        ConShortestPathC(vc,k) = 0;
                        ConShortestPathC(vd,k) = 0;
                    }
                    MS_ShortestPathCollection( va, e1, a );
                    MS_ShortestPathCollection( vd, e2, b );
                    MS_ShortestPathCollection( vc, e4, c );
                    MS_ShortestPathCollection( vb, e5, d );
                    if (ReallyLessCollection(globmin,tmp)) {
                        throw RuntimeException("min increased");
                    }
                    else {
                        globmin = tmp;
                        for(int k=0; k<NT; k++) {
                            globminA[k] = tmpA[k];
                        }
                    }
#endif
                }
            }


        }
    return (nswap);
}

#undef Mab
#undef Mac
#undef Mad
#undef Mae
#undef Mbc
#undef Mbd
#undef Mbe
#undef Mcd
#undef Mce
#undef Mde
#undef Wab
#undef Wac
#undef Wad
#undef Wae
#undef Wbc
#undef Wbd
#undef Wbe
#undef Wcd
#undef Wce
#undef Wde

double MinSquareTreeCollection::FiveSubtreeCollection(double Wab, double Mab, double Wac, double Mac,
                                                      double Wad, double Mad, double Wae, double Mae,
                                                      double Wbc, double Mbc, double Wbd, double Mbd,
                                                      double Wbe, double Mbe, double Wcd, double Mcd,
                                                      double Wce, double Mce, double Wde, double Mde,
                                                      DblVector &L, int Skip) {
    double eps, mmax, T[3];
    int i, it, j;
    DblMatrix Z1(7, 7), Z2(7, 7), Z3(7, 7);
    DblVector B1(7), B2(7), B3(7);

    // incomplete family
    if (Skip > 0) {
        if (Skip == 1) {
            // a missing
            DblVector L4(5);
            double t = FourSubtree(Wbc, Mbc, Wbd, Mbd, Wbe, Mbe, Wcd, Mcd, Wce, Mce, Wde, Mde, L4);
            L[1] = L[2] = L4[0] / 2;
            L[3] = L4[1];
            L[4] = L4[2];
            L[5] = L4[3];
            L[6] = L4[4];
            return (t);
        }
        else if (Skip == 2) {
            // b missing
            DblVector L4(5);
            double t = FourSubtree(Wac, Mac, Wad, Mad, Wae, Mae, Wcd, Mcd, Wce, Mce, Wde, Mde, L4);
            L[0] = L[2] = L4[0] / 2;
            L[3] = L4[1];
            L[4] = L4[2];
            L[5] = L4[3];
            L[6] = L4[4];
            return (t);
        }
        else if (Skip == 4) {
            // c missing
            DblVector L4(5);
            double t = FourSubtree(Wab, Mab, Wad, Mad, Wae, Mae, Wbd, Mbd, Wbe, Mbe, Wde, Mde, L4);
            L[0] = L4[0];
            L[1] = L4[1];
            L[2] = L[4] = L4[2] / 2;
            L[5] = L4[3];
            L[6] = L4[4];
            return (t);
        }
        else if (Skip == 8) {
            // d missing
            DblVector L4(5);
            double t = FourSubtree(Wab, Mab, Wac, Mac, Wae, Mae, Wbc, Mbc, Wbe, Mbe, Wce, Mce, L4);
            L[0] = L4[0];
            L[1] = L4[1];
            L[2] = L4[2];
            L[3] = L4[3];
            L[4] = L[6] = L4[4] / 2;
            return (t);
        }
        else if (Skip == 16) {
            // e missing
            DblVector L4(5);
            double t = FourSubtree(Wab, Mab, Wac, Mac, Wad, Mad, Wbc, Mbc, Wbd, Mbd, Wcd, Mcd, L4);
            L[0] = L4[0];
            L[1] = L4[1];
            L[2] = L4[2];
            L[3] = L4[3];
            L[4] = L[5] = L4[4] / 2;
            return (t);
        }


        else {
            // a,b,c
            if (Skip == 24) {
                ThreeOptimSubset(Wab, Mab, Wac, Mac, Wbc, Mbc, T);
                L[0] = T[0];
                L[1] = T[1];
                L[2] = L[3] = T[2] / 2;
            }
            // a,b,d
            else if (Skip == 20) {
                ThreeOptimSubset(Wab, Mab, Wad, Mad, Wbd, Mbd, T);
                L[0] = T[0];
                L[1] = T[1];
                L[2] = L[4] = L[5] = T[2] / 3;
            }
            // a,b,e
            else if (Skip == 12) {
                ThreeOptimSubset(Wab, Mab, Wae, Mae, Wbe, Mbe, T);
                L[0] = T[0];
                L[1] = T[1];
                L[2] = L[4] = L[6] = T[2] / 3;
            }
            // a,c,d
            else if (Skip == 18) {
                ThreeOptimSubset(Wac, Mac, Wad, Mad, Wcd, Mcd, T);
                L[0] = L[2] = T[0] / 2;
                L[3] = T[1];
                L[4] = L[5] = T[2] / 2;
            }
            // a,c,e
            else if (Skip == 10) {
                ThreeOptimSubset(Wac, Mac, Wae, Mae, Wce, Mce, T);
                L[0] = L[2] = T[0] / 2;
                L[3] = T[1];
                L[4] = L[6] = T[2] / 2;
            }
            // a,d,e
            else if (Skip == 6) {
                ThreeOptimSubset(Wad, Mad, Wae, Mae, Wde, Mde, T);
                L[0] = L[2] = L[4] = T[0] / 3;
                L[5] = T[1];
                L[6] = T[2];
            }
            // b,c,d
            else if (Skip == 17) {
                ThreeOptimSubset(Wbc, Mbc, Wbd, Mbd, Wcd, Mcd, T);
                L[1] = L[2] = T[0] / 2;
                L[3] = T[1];
                L[4] = L[5] = T[2] / 2;
            }
            // b,c,e
            else if (Skip == 9) {
                ThreeOptimSubset(Wbc, Mbc, Wbe, Mbe, Wce, Mce, T);
                L[1] = L[2] = T[0] / 2;
                L[3] = T[1];
                L[4] = L[6] = T[2] / 2;
            }
            // b,d,e
            else if (Skip == 5) {
                ThreeOptimSubset(Wbd, Mbd, Wbe, Mbe, Wde, Mde, T);
                L[1] = L[2] = L[4] = T[0] / 3;
                L[5] = T[1];
                L[6] = T[2];
            }
            // c,d,e
            else if (Skip == 3) {
                ThreeOptimSubset(Wcd, Mcd, Wce, Mce, Wde, Mde, T);
                L[3] = L[4] = T[0] / 2;
                L[5] = T[1];
                L[6] = T[2];
            }
            // a,b
            else if (Skip == 28) {
                L[0] = L[1] = max(0.0, Mab / Wab / 2);
            }
            // a,c
            else if (Skip == 26) {
                L[0] = L[2] = L[3] = max(0.0, Mac / Wac / 3);
            }
            // a,d
            else if (Skip == 22) {
                L[0] = L[2] = L[4] = L[5] = max(0.0, Mad / Wad / 4);
            }
            // a,e
            else if (Skip == 14) {
                L[0] = L[2] = L[4] = L[6] = max(0.0, Mae / Wae / 4);
            }
            // b,c
            else if (Skip == 25) {
                L[1] = L[2] = L[3] = max(0.0, Mbc / Wbc / 3);
            }
            // b,d
            else if (Skip == 21) {
                L[1] = L[2] = L[4] = L[5] = max(0.0, Mbd / Wbd / 4);
            }
            // b,e
            else if (Skip == 13) {
                L[1] = L[2] = L[4] = L[6] = max(0.0, Mbe / Wbe / 4);
            }
            // c,d
            else if (Skip == 19) {
                L[3] = L[4] = L[5] = max(0.0, Mcd / Wcd / 3);
            }
            // c,e
            else if (Skip == 11) {
                L[3] = L[4] = L[6] = max(0.0, Mce / Wce / 3);
            }
            // d,e
            else if (Skip == 7) {
                L[5] = L[6] = max(0.0, Mde / Wde / 2);
            }

            return (0);
        }
    }

    Z1 <<
    Wab + Wac + Wad + Wae,
    Wab,
    Wac + Wad + Wae,
    Wac,
    Wad + Wae,
    Wad,
    Wae,
    Wab,
    Wab + Wbc + Wbd + Wbe,
    Wbc + Wbd + Wbe,
    Wbc,
    Wbd + Wbe,
    Wbd,
    Wbe,
    Wac + Wad + Wae,
    Wbc + Wbd + Wbe,
    Wac + Wad + Wbc + Wbd + Wae + Wbe,
    Wac + Wbc,
    Wad + Wbd + Wae + Wbe,
    Wad + Wbd,
    Wae + Wbe,
    Wac,
    Wbc,
    Wac + Wbc,
    Wac + Wbc + Wcd + Wce,
    Wcd + Wce,
    Wcd,
    Wce,
    Wad + Wae,
    Wbd + Wbe,
    Wad + Wbd + Wae + Wbe,
    Wcd + Wce,
    Wad + Wbd + Wcd + Wae + Wbe + Wce,
    Wad + Wbd + Wcd,
    Wae + Wbe + Wce,
    Wad,
    Wbd,
    Wad + Wbd,
    Wcd,
    Wad + Wbd + Wcd,
    Wad + Wbd + Wcd + Wde,
    Wde,
    Wae,
    Wbe,
    Wae + Wbe,
    Wce,
    Wae + Wbe + Wce,
    Wde,
    Wae + Wbe + Wce + Wde;

    B1 <<
    Mab + Mac + Mad + Mae,
    Mab + Mbc + Mbd + Mbe,
    Mac + Mad + Mbc + Mbd + Mae + Mbe,
    Mac + Mbc + Mcd + Mce,
    Mad + Mbd + Mcd + Mae + Mbe + Mce,
    Mad + Mbd + Mcd + Mde,
    Mae + Mbe + Mce + Mde;

    Z2 = Z1;
    B2 = B1;

    for (it = 0; it < 36; it++) {
        //gausselimi(Z3,B3,L,7);
        // eigen2 syntax:
        //    Z2.lu().solve(B2,&L);
        L = Z2.lu().solve(B2);

        for (mmax = i = 0; i < 7; i++) if (abs(L[i]) > mmax) mmax = abs(L[i]);
        eps = 50 * DBL_EPSILON * mmax;
        for (i = 0; i < 7; i++) {
            if (L[i] < MinLen - eps) {
                double maxcol = {0};
                for (j = 0; j < 7; j++) {
                    if (abs(Z1(j, i)) > maxcol) maxcol = abs(Z1(j, i));
                    Z2(i, j) = 0;
                }
                Z2(i, i) = 10 * maxcol;
                B2[i] = MinLen * 10 * maxcol;
                break;
            } else if (abs(L[i] - MinLen) < eps) {
                double t = -B1[i];
                for (j = 0; j < 7; j++) t += Z1(i, j) * L[j];
                if (t < -eps) {
                    for (j = 0; j < 7; j++) Z2(i, j) = Z1(i, j);
                    B2[i] = B1[i];
                    break;
                }
            }
        }
        if (i >= 7) break;
    }
    for (i = 0; i < 7; i++)
        if (L[i] < MinLen - eps) {
            printf("it=%d, L[%d]=%g, MinLen-L[%d]=%g, eps=%g\n", it, i, L[i],
                   i, MinLen - L[i], eps);
            for (i = 0; i < 49; i++) printf("Z2[%d] := %.18g:\n", i + 1, Z2.data()[i]);
            for (i = 0; i < 7; i++) printf("B2[%d] := %.18g:\n", i + 1, B2.data()[i]);
            printf("MinLen := %.18g:\n", MinLen);
            printf("eps := %.18g:\n", eps);
            printf("mmax := %.18g:\n", mmax);
            throw RuntimeException("internal error in LeastSquaresTree -- length < MinLen in 5-Optim");
        }

    double t1, t4, t7, t10, t13, t16, t19, t22, t29, t36;

    t1 = L[2] + L[4] + L[0] + L[5];
    t4 = L[4] + L[2] + L[1] + L[5];
    t7 = L[2] + L[4] + L[6] + L[0];
    t10 = L[2] + L[4] + L[6] + L[1];
    t13 = L[2] + L[3] + L[0];
    t16 = L[3] + L[2] + L[1];
    t19 = L[4] + L[3] + L[5];
    t22 = L[6] + L[4] + L[3];
    t29 = L[0] + L[1];
    t36 = L[5] + L[6];
    return (t1 * t1 * Wad + t4 * t4 * Wbd + t7 * t7 * Wae + t10 * t10 * Wbe + t13 * t13 * Wac +
            t16 * t16 * Wbc + t19 * t19 * Wcd + t22 * t22 * Wce + t29 * t29 * Wab + t36 * t36 * Wde -
            2 * (t1 * Mad + t4 * Mbd + t7 * Mae + t10 * Mbe + t13 * Mac + t16 * Mbc + t19 * Mcd + t22 * Mce +
                 t29 * Mab + Mde * t36));
}


/***************************************
 #
 #   a              c              d
 #    \             |             /
 #    1\           4|            /6
 #      \     3     |     5     /
 #       o----------o----------o
 #      /x          w          y\
 #    2/                         \7
 #    /                           \
 #   b                             e
 #
 #*************************************/

void MinSquareTreeCollection::SwapFiveSubtree(int e1, int e2, int e4, int e6, int e7,
                                              int ix, int iw, int iy, int e3, int e5) {
    EdgeC[e1].From = ix + ne;
    EdgeC[e2].From = ix + ne;
    EdgeC[e4].From = iw + ne;
    EdgeC[e6].From = iy + ne;
    EdgeC[e7].From = iy + ne;
    EdgeC[e3].From = ix + ne;
    EdgeC[e3].To = iw + ne;
    EdgeC[e5].From = iy + ne;
    EdgeC[e5].To = iw + ne;
    inc(ix, 0) = e1;
    inc(ix, 1) = e2;
    inc(ix, 2) = e3;
    inc(iw, 0) = e3;
    inc(iw, 1) = e4;
    inc(iw, 2) = e5;
    inc(iy, 0) = e5;
    inc(iy, 1) = e6;
    inc(iy, 2) = e7;
}

int MinSquareTreeCollection::DecodeFlags(int f, int *flag) {
    int V[5] = {1, 2, 4, 8, 16};
    flag[0] = (f & V[0]) == V[0];
    flag[1] = (f & V[1]) == V[1];
    flag[2] = (f & V[2]) == V[2];
    flag[3] = (f & V[3]) == V[3];
    flag[4] = (f & V[4]) == V[4];
    return (flag[0] + flag[1] + flag[2] + flag[3] + flag[4]);
}

int MinSquareTreeCollection::FiveOptimCollectionExtended345(int e3, int e4, int e5, int vc, int ix,
                                                            int iy, int iw, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W, int *ptrnswap) {
    double t1, t2;
    int e1, e2, e6, e7, i, va, vb, vd, ve, flag[8], tmp, done;
    int NewTopology, E[5];
    const int C[15][5] = {
        {0, 1, 2, 3, 4},
        {0, 3, 2, 1, 4},
        {0, 4, 2, 3, 1},
        {0, 2, 1, 3, 4},
        {0, 3, 1, 2, 4},
        {0, 4, 1, 2, 3},
        {0, 1, 3, 2, 4},
        {0, 2, 3, 1, 4},
        {0, 4, 3, 1, 2},
        {0, 1, 4, 2, 3},
        {0, 2, 4, 1, 3},
        {0, 3, 4, 1, 2},
        {1, 2, 0, 3, 4},
        {2, 3, 0, 1, 4},
        {2, 4, 0, 1, 3},
    };

    if (inc(ix, 0) == e3) {
        e1 = inc(ix, 1);
        e2 = inc(ix, 2);
    }
    else if (inc(ix, 1) == e3) {
        e1 = inc(ix, 0);
        e2 = inc(ix, 2);
    }
    else {
        e1 = inc(ix, 0);
        e2 = inc(ix, 1);
    }

    if (inc(iy, 0) == e5) {
        e6 = inc(iy, 1);
        e7 = inc(iy, 2);
    }
    else if (inc(iy, 1) == e5) {
        e6 = inc(iy, 0);
        e7 = inc(iy, 2);
    }
    else {
        e6 = inc(iy, 0);
        e7 = inc(iy, 1);
    }

    /*e1*/
    if (EdgeC[e1].From == ix + ne) {
        va = EdgeC[e1].To;
    }
    else {
        va = EdgeC[e1].From;
        EdgeC[e1].From = EdgeC[e1].To;
        EdgeC[e1].To = va;
    }
    /*e2*/
    if (EdgeC[e2].From == ix + ne) {
        vb = EdgeC[e2].To;
    }
    else {
        vb = EdgeC[e2].From;
        EdgeC[e2].From = EdgeC[e2].To;
        EdgeC[e2].To = vb;
    }
    /*e4*/
    if (EdgeC[e4].From == vc) {
        EdgeC[e4].From = EdgeC[e4].To;
        EdgeC[e4].To = vc;
    }
    /*e6*/
    if (EdgeC[e6].From == iy + ne) {
        vd = EdgeC[e6].To;
    }
    else {
        vd = EdgeC[e6].From;
        EdgeC[e6].From = EdgeC[e6].To;
        EdgeC[e6].To = vd;
    }
    /*e7*/
    if (EdgeC[e7].From == iy + ne) {
        ve = EdgeC[e7].To;
    }
    else {
        ve = EdgeC[e7].From;
        EdgeC[e7].From = EdgeC[e7].To;
        EdgeC[e7].To = ve;
    }

    for (int k = 0; k < NT; k++) {
        ConShortestPathC(va, k) = 0;
        ConShortestPathC(vb, k) = 0;
        ConShortestPathC(vc, k) = 0;
        ConShortestPathC(vd, k) = 0;
        ConShortestPathC(ve, k) = 0;
    }
    MS_ShortestPathCollection(va, e1, a);
    MS_ShortestPathCollection(vb, e2, b);
    MS_ShortestPathCollection(vc, e4, c);
    MS_ShortestPathCollection(vd, e6, d);
    MS_ShortestPathCollection(ve, e7, e);

    /* in Skip integer: a->1, b->2, c->4, d->8, e->16 */
    for (int k = 0; k < NT; k++) {

        flag[a] = flag[b] = flag[c] = flag[d] = flag[e] = 0;

        for (int j = 0; j < ne; j++)
            if (!isNA(k, j))
                flag[ShortestLabel[j]] = 1;

        Skip[k] = 0;
        Skip[k] += flag[a] ? 0 : 1;
        Skip[k] += flag[b] ? 0 : 2;
        Skip[k] += flag[c] ? 0 : 4;
        Skip[k] += flag[d] ? 0 : 8;
        Skip[k] += flag[e] ? 0 : 16;
    }

    E[0] = e1;
    E[1] = e2;
    E[2] = e4;
    E[3] = e6;
    E[4] = e7;

    t1 = DistanceFitCollection();
    NewTopology = 0;
    done = 0;
    i = 0;
    while (1) {
        SwapFiveSubtree(E[C[i][0]], E[C[i][1]], E[C[i][2]], E[C[i][3]], E[C[i][4]],
                        ix, iw, iy, e3, e5);
        /* now select the edges to be fitted */
        for (int k = 0; k < NT; k++) {
            /* ok, now we need to find out which edge need to be fitted */
            tmp = DecodeFlags(Skip[k], flag);
            if (tmp <= 3) {
                for (int j = 0; j < 5; j++) {
                    if (flag[C[i][j]] == 0) {
                        EdgeC[E[C[i][j]]].len[k] = -EdgeC[E[C[i][j]]].len[k];
                    }
                }
                if (flag[C[i][0]] == 0 || flag[C[i][1]] == 0) {
                    if (EdgeC[e3].len[k] == DBL_MAX) {
                        EdgeC[e3].len[k] = MinLen;
                    }
                    EdgeC[e3].len[k] = -EdgeC[e3].len[k];
                }
                if (flag[C[i][3]] == 0 || flag[C[i][4]] == 0) {
                    /* we use the consensus length as "starting point" */
                    if (EdgeC[e5].len[k] == DBL_MAX) {
                        EdgeC[e5].len[k] = MinLen;
                    }
                    EdgeC[e5].len[k] = -EdgeC[e5].len[k];
                }
            }
        }

        FitLabeledEdgesC(0);

        if (done >= 1) {
            t2 = DistanceFitCollection();
            if (!ReallyLessCollection(t1, t2)) {
                printf("done; t=%g\n", t2);
                break;
            }
            else {
                printf("Still suboptimal=%d,NewTop=%d,t1=%g,t2=%g\n", i,
                       NewTopology, t1, t2);
                if (done > 3) {
                    printf("Giving up: we reoptimize all edges\n");
                    FitLabeledEdgesC(1);
                }
                done++;
                continue;
            }
        }

        t2 = DistanceFitCollection();

        if (ReallyLessCollection(t2, t1)) {
            t1 = t2;
            NewTopology = i;
        }
        i++;
        if (i == 15) {
            done = 1;
            i = NewTopology;
        }
    }

    if (i > 0)
        (*ptrnswap)++;
    return (NewTopology);
}

inline int MinSquareTreeCollection::ModifySkip(int flag, int p0, int p1, int p2, int p3, int p4) {
    int V[5] = {1, 2, 4, 8, 16};
    int res, t[5];
    //decode
    if (flag == 0)
        return (0);
    t[0] = (flag & V[0]) == V[0];
    t[1] = (flag & V[1]) == V[1];
    t[2] = (flag & V[2]) == V[2];
    t[3] = (flag & V[3]) == V[3];
    t[4] = (flag & V[4]) == V[4];

    res = V[0] * t[p0]
    + V[1] * t[p1]
    + V[2] * t[p2]
    + V[3] * t[p3]
    + V[4] * t[p4];
    return (res);
}

/* the following sets EdgeC[k].len[eX] to 0 for all eX until a split
 with two existing branches is encountered

 /
 (missing, i.e. len=DBL_MAX)
 -> start      /
 ________.
 \
 \           /
 \         /
 '-------*   <- end
 |        \
 (missing)           \
 |


 */
void MinSquareTreeCollection::delPathLength(int k, int from, int ExcludedEdge) {
    int dest, tmp, t;

    if (from < ne)
        return;

    tmp = 0;
    t = -1;
    for (int j = 0; j < 3; j++) {
        int aE = inc(from - ne, j);
        if (aE != ExcludedEdge && EdgeC[aE].len[k] != DBL_MAX) {
            tmp++;
            t = aE;
        }
    }
    if (tmp == 2 || tmp == 0) {
        return;
    }
    else {
        EdgeC[t].len[k] = 0;
        if (EdgeC[t].From == from) {
            dest = EdgeC[t].To;
        }
        else {
            dest = EdgeC[t].From;
        }
        delPathLength(k, dest, t);
    }

}

int MinSquareTreeCollection::FiveOptimCollection345(int DoNotModify, int e3, int e4, int e5, int vc, int ix,
                                                    int iy, int iw, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W, int *ptrnswap)
/* requires inc[][] vector */
{
    double dd, t1, t2, Wei;
    int e1, e2, e6, e7, va, vb, vd, ve, flag[8], tot;
    int NewTopology;
    DblVector Lout = DblVector::Zero(7);

#ifdef DEBUG
    double tmp;
#endif

    if (inc(ix, 0) == e3) {
        e1 = inc(ix, 1);
        e2 = inc(ix, 2);
    }
    else if (inc(ix, 1) == e3) {
        e1 = inc(ix, 0);
        e2 = inc(ix, 2);
    }
    else {
        e1 = inc(ix, 0);
        e2 = inc(ix, 1);
    }

    if (inc(iy, 0) == e5) {
        e6 = inc(iy, 1);
        e7 = inc(iy, 2);
    }
    else if (inc(iy, 1) == e5) {
        e6 = inc(iy, 0);
        e7 = inc(iy, 2);
    }
    else {
        e6 = inc(iy, 0);
        e7 = inc(iy, 1);
    }

    va = EdgeC[e1].From == ix + ne ? EdgeC[e1].To : EdgeC[e1].From;
    vb = EdgeC[e2].From == ix + ne ? EdgeC[e2].To : EdgeC[e2].From;
    vd = EdgeC[e6].From == iy + ne ? EdgeC[e6].To : EdgeC[e6].From;
    ve = EdgeC[e7].From == iy + ne ? EdgeC[e7].To : EdgeC[e7].From;

    LabelNonExistEdges();

    for (int k = 0; k < NT; k++) {
        ConShortestPathC(va, k) = 0;
        ConShortestPathC(vb, k) = 0;
        ConShortestPathC(vc, k) = 0;
        ConShortestPathC(vd, k) = 0;
        ConShortestPathC(ve, k) = 0;
        delPathLength(k, va, e1);
        delPathLength(k, vb, e2);
        delPathLength(k, vc, e4);
        delPathLength(k, vd, e6);
        delPathLength(k, ve, e7);

    }

    MS_ShortestPathCollection(va, e1, a);
    MS_ShortestPathCollection(vb, e2, b);
    MS_ShortestPathCollection(vc, e4, c);
    MS_ShortestPathCollection(vd, e6, d);
    MS_ShortestPathCollection(ve, e7, e);


    /* in Skip integer: a=1, b=2, c=4, d=8, e=16 */
    for (int k = 0; k < NT; k++) {
        for (int i = 0; i <= 11; i++)
            W(k, i) = M(k, i) = 0.0;

        flag[a] = flag[b] = flag[c] = flag[d] = flag[e] = 0;

        for (int si = 0; si < ne; si++) {
            if (isNA(k, si))
                continue;
            flag[ShortestLabel[si]]++;

            for (int sj = si + 1; sj < ne; sj++) {
                if (ShortestLabel[si] != ShortestLabel[sj]) {
                    if (isNA(k, sj))
                        continue;
                    int i = ShortestLabel[si] + ShortestLabel[sj];
                    double v = aVAR(k, si, sj);
                    if (v == 0) {
                        Wei = 1.0e-10;
                        dd = 300.0;
                    }
                    else {
                        Wei = v;
                        dd = aDIST(k, si, sj);
                    }
                    W(k, i) += Wei;
                    M(k, i) += Wei * (dd - ConShortestPathC(si, k)
                                      - ConShortestPathC(sj, k));
                }
            }
        }
        Skip[k] = 0;
        Skip[k] += flag[a] ? 0 : 1;
        Skip[k] += flag[b] ? 0 : 2;
        Skip[k] += flag[c] ? 0 : 4;
        Skip[k] += flag[d] ? 0 : 8;
        Skip[k] += flag[e] ? 0 : 16;
    }
    tot = 0;
    for (int k = 0; k < NT; k++)
        if (Skip[k] == 0)
            tot++;
    if (tot > 0) {
        EdgeC[e1].len0 = 0;
        EdgeC[e2].len0 = 0;
        EdgeC[e3].len0 = 0;
        EdgeC[e4].len0 = 0;
        EdgeC[e5].len0 = 0;
        EdgeC[e6].len0 = 0;
        EdgeC[e7].len0 = 0;
    }

    NewTopology = 0;
    t1 = 0;
    for (int k = 0; k < NT; k++) {
        t1 += FiveSubtreeCollection(W(k, a + b), M(k, a + b), W(k, a + c), M(k, a + c),
                                    W(k, a + d), M(k, a + d), W(k, a + e), M(k, a + e), W(k, b + c), M(k, b + c),
                                    W(k, b + d), M(k, b + d), W(k, b + e), M(k, b + e), W(k, c + d), M(k, c + d),
                                    W(k, c + e), M(k, c + e), W(k, d + e), M(k, d + e), Lout, Skip[k]);
        L.row(k) = Lout;
        EdgeC[e1].len[k] = L(k, 0);
        EdgeC[e2].len[k] = L(k, 1);
        EdgeC[e3].len[k] = L(k, 2);
        EdgeC[e4].len[k] = L(k, 3);
        EdgeC[e5].len[k] = L(k, 4);
        EdgeC[e6].len[k] = L(k, 5);
        EdgeC[e7].len[k] = L(k, 6);

        /* replacement idea:
         */

        if (Skip[k] == 0) {
            EdgeC[e1].len0 += L(k, 0) / tot;
            EdgeC[e2].len0 += L(k, 1) / tot;
            EdgeC[e3].len0 += L(k, 2) / tot;
            EdgeC[e4].len0 += L(k, 3) / tot;
            EdgeC[e5].len0 += L(k, 4) / tot;
            EdgeC[e6].len0 += L(k, 5) / tot;
            EdgeC[e7].len0 += L(k, 6) / tot;
        }
    }

    if (DoNotModify)
        return (0);

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + d), M(k, a + d), W(k, a + c), M(k, a + c),
                                    W(k, a + b), M(k, a + b), W(k, a + e), M(k, a + e), W(k, c + d), M(k, c + d),
                                    W(k, b + d), M(k, b + d), W(k, d + e), M(k, d + e), W(k, b + c), M(k, b + c),
                                    W(k, c + e), M(k, c + e), W(k, b + e), M(k, b + e), Lout,
                                    ModifySkip(Skip[k], 0, 3, 2, 1, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 += L(k, 0) / tot;
                EdgeC[e2].len0 += L(k, 1) / tot;
                EdgeC[e3].len0 += L(k, 2) / tot;
                EdgeC[e4].len0 += L(k, 3) / tot;
                EdgeC[e5].len0 += L(k, 4) / tot;
                EdgeC[e6].len0 += L(k, 5) / tot;
                EdgeC[e7].len0 += L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vd;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vc;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vb;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 1;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 1, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif

    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + e), M(k, a + e), W(k, a + c), M(k, a + c), W(k, a + d), M(k, a + d), W(k, a + b), M(k, a + b), W(k, c + e), M(k, c + e), W(k, d + e), M(k, d + e), W(k, b + e), M(k, b + e), W(k, c + d), M(k, c + d), W(k, b + c), M(k, b + c), W(k, b + d), M(k, b + d), Lout, ModifySkip(Skip[k], 0, 4, 2, 3, 1));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 += L(k, 0) / tot;
                EdgeC[e2].len0 += L(k, 1) / tot;
                EdgeC[e3].len0 += L(k, 2) / tot;
                EdgeC[e4].len0 += L(k, 3) / tot;
                EdgeC[e5].len0 += L(k, 4) / tot;
                EdgeC[e6].len0 += L(k, 5) / tot;
                EdgeC[e7].len0 += L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = ve;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vc;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vb;
        NewTopology = 2;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 2, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, b + c), M(k, b + c), W(k, a + c), M(k, a + c),
                                    W(k, c + d), M(k, c + d), W(k, c + e), M(k, c + e),
                                    W(k, a + b), M(k, a + b), W(k, b + d), M(k, b + d), W(k, b + e), M(k, b + e),
                                    W(k, a + d), M(k, a + d), W(k, a + e), M(k, a + e), W(k, d + e), M(k, d + e), Lout,
                                    ModifySkip(Skip[k], 2, 1, 0, 3, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 += L(k, 0) / tot;
                EdgeC[e2].len0 += L(k, 1) / tot;
                EdgeC[e3].len0 += L(k, 2) / tot;
                EdgeC[e4].len0 += L(k, 3) / tot;
                EdgeC[e5].len0 += L(k, 4) / tot;
                EdgeC[e6].len0 += L(k, 5) / tot;
                EdgeC[e7].len0 += L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = vc;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vb;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = va;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 3;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 3, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, c + d), M(k, c + d), W(k, a + c), M(k, a + c), W(k, b + c), M(k, b + c), W(k, c + e), M(k, c + e), W(k, a + d), M(k, a + d), W(k, b + d), M(k, b + d), W(k, d + e), M(k, d + e), W(k, a + b), M(k, a + b), W(k, a + e), M(k, a + e), W(k, b + e), M(k, b + e), Lout, ModifySkip(Skip[k], 2, 3, 0, 1, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 += L(k, 0) / tot;
                EdgeC[e2].len0 += L(k, 1) / tot;
                EdgeC[e3].len0 += L(k, 2) / tot;
                EdgeC[e4].len0 += L(k, 3) / tot;
                EdgeC[e5].len0 += L(k, 4) / tot;
                EdgeC[e6].len0 += L(k, 5) / tot;
                EdgeC[e7].len0 += L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = vc;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vd;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = va;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vb;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 4;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 4, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, c + e), M(k, c + e), W(k, a + c), M(k, a + c),
                                    W(k, c + d), M(k, c + d), W(k, b + c), M(k, b + c),
                                    W(k, a + e), M(k, a + e), W(k, d + e), M(k, d + e), W(k, b + e), M(k, b + e),
                                    W(k, a + d), M(k, a + d), W(k, a + b), M(k, a + b), W(k, b + d), M(k, b + d), Lout,
                                    ModifySkip(Skip[k], 2, 4, 0, 3, 1));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 += L(k, 0) / tot;
                EdgeC[e2].len0 += L(k, 1) / tot;
                EdgeC[e3].len0 += L(k, 2) / tot;
                EdgeC[e4].len0 += L(k, 3) / tot;
                EdgeC[e5].len0 += L(k, 4) / tot;
                EdgeC[e6].len0 += L(k, 5) / tot;
                EdgeC[e7].len0 += L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = vc;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = ve;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = va;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vb;
        NewTopology = 5;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 5, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + c), M(k, a + c), W(k, a + b), M(k, a + b),
                                    W(k, a + d), M(k, a + d), W(k, a + e), M(k, a + e),
                                    W(k, b + c), M(k, b + c), W(k, c + d), M(k, c + d), W(k, c + e), M(k, c + e),
                                    W(k, b + d), M(k, b + d), W(k, b + e), M(k, b + e), W(k, d + e), M(k, d + e), Lout,
                                    ModifySkip(Skip[k], 0, 2, 1, 3, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vc;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vb;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 6;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 6, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + d), M(k, a + d), W(k, a + b), M(k, a + b),
                                    W(k, a + c), M(k, a + c), W(k, a + e), M(k, a + e),
                                    W(k, b + d), M(k, b + d), W(k, c + d), M(k, c + d), W(k, d + e), M(k, d + e),
                                    W(k, b + c), M(k, b + c), W(k, b + e), M(k, b + e), W(k, c + e), M(k, c + e), Lout,
                                    ModifySkip(Skip[k], 0, 3, 1, 2, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vd;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vb;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vc;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 7;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 7, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + e), M(k, a + e), W(k, a + b), M(k, a + b),
                                    W(k, a + d), M(k, a + d), W(k, a + c), M(k, a + c),
                                    W(k, b + e), M(k, b + e), W(k, d + e), M(k, d + e), W(k, c + e), M(k, c + e),
                                    W(k, b + d), M(k, b + d), W(k, b + c), M(k, b + c), W(k, c + d), M(k, c + d), Lout,
                                    ModifySkip(Skip[k], 0, 4, 1, 3, 2));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = ve;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vb;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vc;
        NewTopology = 8;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 8, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + b), M(k, a + b), W(k, a + d), M(k, a + d),
                                    W(k, a + c), M(k, a + c), W(k, a + e), M(k, a + e),
                                    W(k, b + d), M(k, b + d), W(k, b + c), M(k, b + c), W(k, b + e), M(k, b + e),
                                    W(k, c + d), M(k, c + d), W(k, d + e), M(k, d + e), W(k, c + e), M(k, c + e), Lout,
                                    ModifySkip(Skip[k], 0, 1, 3, 2, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vb;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vd;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vc;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 9;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 9, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + c), M(k, a + c), W(k, a + d), M(k, a + d),
                                    W(k, a + b), M(k, a + b), W(k, a + e), M(k, a + e),
                                    W(k, c + d), M(k, c + d), W(k, b + c), M(k, b + c), W(k, c + e), M(k, c + e),
                                    W(k, b + d), M(k, b + d), W(k, d + e), M(k, d + e), W(k, b + e), M(k, b + e), Lout,
                                    ModifySkip(Skip[k], 0, 2, 3, 1, 4));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vc;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vd;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vb;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = ve;
        NewTopology = 10;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 10, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + e), M(k, a + e), W(k, a + d), M(k, a + d),
                                    W(k, a + c), M(k, a + c), W(k, a + b), M(k, a + b),
                                    W(k, d + e), M(k, d + e), W(k, c + e), M(k, c + e), W(k, b + e), M(k, b + e),
                                    W(k, c + d), M(k, c + d), W(k, b + d), M(k, b + d), W(k, b + c), M(k, b + c), Lout,
                                    ModifySkip(Skip[k], 0, 4, 3, 2, 1));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = ve;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = vd;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vc;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vb;
        NewTopology = 11;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 11, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + b), M(k, a + b), W(k, a + e), M(k, a + e),
                                    W(k, a + d), M(k, a + d), W(k, a + c), M(k, a + c),
                                    W(k, b + e), M(k, b + e), W(k, b + d), M(k, b + d), W(k, b + c), M(k, b + c),
                                    W(k, d + e), M(k, d + e), W(k, c + e), M(k, c + e), W(k, c + d), M(k, c + d), Lout,
                                    ModifySkip(Skip[k], 0, 1, 4, 3, 2));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vb;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = ve;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vc;
        NewTopology = 12;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 12, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + d), M(k, a + d), W(k, a + e), M(k, a + e),
                                    W(k, a + b), M(k, a + b), W(k, a + c), M(k, a + c),
                                    W(k, d + e), M(k, d + e), W(k, b + d), M(k, b + d), W(k, c + d), M(k, c + d),
                                    W(k, b + e), M(k, b + e), W(k, c + e), M(k, c + e), W(k, b + c), M(k, b + c), Lout,
                                    ModifySkip(Skip[k], 0, 3, 4, 1, 2));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vd;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = ve;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vb;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vc;
        NewTopology = 13;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 13, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }

    t2 = 0;
    for (int k = 0; k < NT; k++) {
        t2 += FiveSubtreeCollection(W(k, a + c), M(k, a + c), W(k, a + e), M(k, a + e),
                                    W(k, a + d), M(k, a + d), W(k, a + b), M(k, a + b),
                                    W(k, c + e), M(k, c + e), W(k, c + d), M(k, c + d), W(k, b + c), M(k, b + c),
                                    W(k, d + e), M(k, d + e), W(k, b + e), M(k, b + e), W(k, b + d), M(k, b + d), Lout,
                                    ModifySkip(Skip[k], 0, 2, 4, 3, 1));
        L.row(k) = Lout;
    }
    if (ReallyLessCollection(t2, t1)) {
        t1 = t2;
        (*ptrnswap)++;
        for (int k = 0; k < NT; k++) {
            EdgeC[e1].len[k] = L(k, 0);
            EdgeC[e2].len[k] = L(k, 1);
            EdgeC[e3].len[k] = L(k, 2);
            EdgeC[e4].len[k] = L(k, 3);
            EdgeC[e5].len[k] = L(k, 4);
            EdgeC[e6].len[k] = L(k, 5);
            EdgeC[e7].len[k] = L(k, 6);
            if (Skip[k] == 0) {
                EdgeC[e1].len0 = L(k, 0) / tot;
                EdgeC[e2].len0 = L(k, 1) / tot;
                EdgeC[e3].len0 = L(k, 2) / tot;
                EdgeC[e4].len0 = L(k, 3) / tot;
                EdgeC[e5].len0 = L(k, 4) / tot;
                EdgeC[e6].len0 = L(k, 5) / tot;
                EdgeC[e7].len0 = L(k, 6) / tot;
            }
        }
        EdgeC[e1].From = ix + ne;
        EdgeC[e1].To = va;
        EdgeC[e2].From = ix + ne;
        EdgeC[e2].To = vc;
        EdgeC[e4].From = iw + ne;
        EdgeC[e4].To = ve;
        EdgeC[e6].From = iy + ne;
        EdgeC[e6].To = vd;
        EdgeC[e7].From = iy + ne;
        EdgeC[e7].To = vb;
        NewTopology = 14;
#ifdef DEBUG
        IncidencesC();
        tmp = DistanceFitCollection();
        printf("FiveOptim: 14, edge %d %.3f\n",e3,tmp);

        if (ReallyLessCollection(globmin,tmp)) {
            printf("  warning: min increased from %.3f\n",globmin);
        }
        else {
            globmin = tmp;
            for(int k=0; k<NT; k++) {
                globminA[k] = tmpA[k];
            }
        }
#endif
    }
    return (NewTopology);
}

int MinSquareTreeCollection::FiveOptimCollection(int mode, DblMatrix &L, IntVector &Skip, DblMatrix &M, DblMatrix &W) {
    int b1, e3, e4, e5, i, ix, iy, vc, nswap, DoNotModify;
    IntVector DoEdge = IntVector::Zero(ne - 2);

    /***************************************
     #
     #   a              c              d
     #    \             |             /
     #    1\           4|            /6
     #      \     3     |     5     /
     #       o----------o----------o
     #      /x          w          y\
     #    2/                         \7
     #    /                           \
     #   b                             e
     #
     #*************************************/
    nswap = 0;

    DoNotModify = mode == 1 ? 1 : 0;

    IncidencesC();
#ifdef DEBUG
    double d1 = DistanceFitCollection();
    if (ReallyLessCollection(globmin,d1)) {
        printf("min increased: from %.3f to %.3f\n",globmin,d1);
    }
    else {
        globmin = d1;
        for(int k=0; k<NT; k++) {
            globminA[k] = tmpA[k];
        }
    }
#endif


    for (i = ne - 3; i >= 0; i--) DoEdge[i] = 1;

    /* in very rare ill-conditioned cases, it may loop despite of
     * ReallyLessCollection(a,b).  Hence we limit the number of iterations
     * over the whole set of internal nodes to 5 */
    for (int iter = 0; iter < 5; iter++) {
        for (i = ne - 3; i >= 0 && !DoEdge[i]; i--);
        if (i < 0) break;  /* no more edges to process */

        for (int iw = ne - 3; iw >= 0; iw--)
            if (DoEdge[iw]) {
                DoEdge[iw] = 0;
                e3 = inc(iw, 0);
                e4 = inc(iw, 1);
                e5 = inc(iw, 2);

                if (EdgeC[e3].From < ne || EdgeC[e3].To < ne) {
                    if (EdgeC[e4].From >= ne && EdgeC[e4].To >= ne) {
                        i = e4;
                        e4 = e3;
                        e3 = i;
                    }
                    else continue;
                }
                if (EdgeC[e5].From < ne || EdgeC[e5].To < ne) {
                    if (EdgeC[e4].From >= ne && EdgeC[e4].To >= ne) {
                        i = e4;
                        e4 = e5;
                        e5 = i;
                    }
                    else continue;
                }

                vc = EdgeC[e4].From == iw + ne ? EdgeC[e4].To : EdgeC[e4].From;
                ix = (EdgeC[e3].From == iw + ne ? EdgeC[e3].To : EdgeC[e3].From) - ne;
                iy = (EdgeC[e5].From == iw + ne ? EdgeC[e5].To : EdgeC[e5].From) - ne;

                if (mode == 2)
                    b1 = FiveOptimCollectionExtended345(e3, e4, e5, vc, ix, iy, iw,
                                                        L, Skip, M, W, &nswap);
                else
                    b1 = FiveOptimCollection345(DoNotModify, e3, e4, e5, vc, ix, iy, iw,
                                                L, Skip, M, W, &nswap);

                if (!b1 && vc >= ne) {
                    if (mode == 2)
                        b1 = FiveOptimCollectionExtended345(e4, e3, e5, ix + ne, vc - ne,
                                                            iy, iw, L, Skip, M, W, &nswap);
                    else
                        b1 = FiveOptimCollection345(DoNotModify, e4, e3, e5, ix + ne, vc - ne,
                                                    iy, iw, L, Skip, M, W, &nswap);
                    if (b1) b1 += 100;
                }
                if (!b1 && vc >= ne) {
                    if (mode == 2)
                        b1 = FiveOptimCollectionExtended345(e3, e5, e4, iy + ne, ix,
                                                            vc - ne, iw, L, Skip, M, W, &nswap);
                    else
                        b1 = FiveOptimCollection345(DoNotModify, e3, e5, e4, iy + ne, ix,
                                                    vc - ne, iw, L, Skip, M, W, &nswap);
                    if (b1) b1 += 200;
                }

                if (b1) {

                    IncidencesC();
                    DoEdge[iw] = DoEdge[ix] = DoEdge[iy] = 1;
#ifdef DEBUG
                    double d2 = DistanceFitCollection();
                    if (ReallyLessCollection(globmin,d2)) {
                        printf("min increased: from %.3f to %.3f\n",globmin,d2);
                    }
                    else {
                        globmin = d2;
                        for(int k=0; k<NT; k++) {
                            globminA[k] = tmpA[k];
                        }
                    }
                    if( ReallyLessCollection(d1,d2) )
                        printf( "LeastSquaresTree increases error - initial tree may not respect MinLen\n");
                    printf("FiveOptim: %d, node %d, lse=%.15g\n", b1,iw,d2);

#endif
                }
            }
    }

    return (nswap);
}

int MinSquareTreeCollection::MapTree2InternalC(const PhyTree &t) throw(ParameterException) {
    int e1, e2, in1, in2, in3;
    if (t.isLeaf()) {
        for (int i = 0; i < ne; i++) {
            if (t.getName() == Labels[i])
                return (i);
        }
        throw ParameterException("label in Initial tree not found in label list: " + t.getName());
    }
    if (t.n_children() != 2) {
        throw ParameterException("Tree node is too short");
    }

    in3 = NewInternalNode++;
    in1 = MapTree2InternalC(t[0]);
    in2 = MapTree2InternalC(t[1]);

    inc(in3, 0) = e1 = NewEdgeIndex++;
    EdgeC[e1].From = in1;
    EdgeC[e1].To = in3 + ne;
    EdgeC[e1].len0 = t[0].getBranchLength();
    for (int k = 0; k < NT; k++) {
        EdgeC[e1].len[k] = t[0].getBranchLength();
    }
    if (in1 >= ne) inc(in1 - ne, 2) = e1;

    inc(in3, 1) = e2 = NewEdgeIndex++;
    EdgeC[e2].From = in2;
    EdgeC[e2].To = in3 + ne;
    EdgeC[e2].len0 = t[1].getBranchLength();
    for (int k = 0; k < NT; k++) {
        EdgeC[e2].len[k] = t[1].getBranchLength();
    }
    if (in2 >= ne) inc(in2 - ne, 2) = e2;

    return (in3 + ne);
}

double MinSquareTreeCollection::EdgeCLength(int i) {
    int ct;
    double t;

    return (EdgeC[i].len0);

    /**** JUNK BELOW
     t = ct = 0;
     for (int k = 0; k < NT; k++) {
     if (EdgeC[i].len[k] != EdgeC[i].len0 && EdgeC[i].len[k] > 0) {
     t += EdgeC[i].len[k];
     ct++;
     }
     }
     if (ct >= 1) return (t / ct);
     else return (EdgeC[i].len0);
     ****/
}

template<class T>
static string to_str(T x) {
    stringstream ss;
    ss << x;
    return ss.str();
}

PhyTree::TreePtr MinSquareTreeCollection::MST_TreeCR(int node, int edgefrom) {
    int e1, e2;

    if (node < ne)
        return (make_shared<PhyTree>(Labels.size() == 0 ? to_str(node + 1) : Labels[node]));

    if (inc(node - ne, 0) == edgefrom) {
        e1 = inc(node - ne, 1);
        e2 = inc(node - ne, 2);
    }
    else if (inc(node - ne, 1) == edgefrom) {
        e1 = inc(node - ne, 0);
        e2 = inc(node - ne, 2);
    }
    else {
        e1 = inc(node - ne, 0);
        e2 = inc(node - ne, 1);
    }


    auto tree = make_shared<PhyTree>();
    tree->addChild(MST_TreeCR(EdgeC[e1].From == node ? EdgeC[e1].To : EdgeC[e1].From, e1), EdgeCLength(e1));
    tree->addChild(MST_TreeCR(EdgeC[e2].From == node ? EdgeC[e2].To : EdgeC[e2].From, e2), EdgeCLength(e2));

    return (tree);

}

/* Build a Darwin Tree from Edge list */
PhyTree::TreePtr MinSquareTreeCollection::getPhyTree() {
    int e1, n[2], bestE = 0;
    double p[2], bestDiff;

    IncidencesC();
    bestDiff = DBL_MAX;
    for (e1 = 0; e1 < 2 * ne - 3; e1++) {
        for (int k = 0; k < NT; k++) {
            ConShortestPathC(EdgeC[e1].To, k) = 0;
            ConShortestPathC(EdgeC[e1].From, k) = 0;
        }
        MS_ShortestPathCollection(EdgeC[e1].To, e1, 0);
        MS_ShortestPathCollection(EdgeC[e1].From, e1, 1);
        p[0] = p[1] = 0;
        n[0] = n[1] = 0;
        /* WARNING: we root the tree based on first tree (family) only */
        for (int s0 = 0; s0 < ne; s0++) {
            n[ShortestLabel[s0]]++;
            p[ShortestLabel[s0]] += ConShortestPathC(s0, 0);
        }
        if (!(n[0] > 0 && n[1] > 0)) throw RuntimeException("assertion failed: n[0] > 0 && n[1] > 0");
        if (abs(p[0] / n[0] - p[1] / n[1]) < bestDiff) {
            bestE = e1;
            bestDiff = abs(p[0] / n[0] - p[1] / n[1]);
        }
        /* if( abs(p[0]/n[0]-p[1]/n[1]) <= EdgeCLength(e1) ) break; */
    }
    if (!(bestE < 2 * ne - 3)) throw RuntimeException("assertion failed: bestE < 2*ne-3");
    e1 = bestE;
    PhyTree::TreePtr tree = make_shared<PhyTree>();
    tree->addChild(MST_TreeCR(EdgeC[e1].From, e1), EdgeCLength(e1) / 2.0);
    tree->addChild(MST_TreeCR(EdgeC[e1].To, e1), EdgeCLength(e1) / 2.0);

    return (tree);
}


MinSquareTreeCollection::MinSquareTreeCollection(const vector<DblMatrix> &matrices, const IntMatrix &mapping, const vector<string> &labels, const PhyTree &tree) throw(ParameterException) {
    initialized = false;
    computed = false;

    //    Edge = BestEdge = NULL;
    //    EdgeC = NULL;

    aDistVar = matrices;
    NT = aDistVar.size();
    aMap = mapping;
    if (NT != aMap.rows()) {
        throw ParameterException("the arrays of distance and vars have different lengths");
    }


    ne = aMap.cols();

    if ((int) labels.size() != 0 && (int) labels.size() != ne) {
        throw ParameterException("incorrect dimensions in arguments for MinSquareTree");
    }
    Labels = labels;

    for (int k = 0; k < NT; k++) {
        /* verify the bounds of mapping array */
        int mmax = -1;
        for (int j = 0; j < ne; j++) {
            if (aMap(k, j) > mmax)
                mmax = aMap(k, j);
            if (aMap(k, j) == 0 || aMap(k, j) > ne || aMap(k, j) < -1) {
                throw ParameterException("mapping array has value outside range");
            }
        }
        if (mmax > ne || aDistVar[k].rows() != aDistVar[k].cols() || aDistVar[k].rows() != mmax) {
            throw ParameterException("mapping array has inconsistent entries");
        }
    }

    MinLen = 1e-6; // Reduced MinLen from 1e-2 - KG

    /* allocate working storage */
    /* EdgeC.len for all Edges */

    /* verify all variances */
    for (int k = 0; k < NT; k++) {
        if ((aDistVar[k].array() < 0).any()) {
            throw ParameterException("distances/variances cannot be < 0");
        }
    }

    /* process additional arguments */
    NewInternalNode = NewEdgeIndex = 0;
    for (int i = 0; i < ne; i++)
        for (int j = i + 1; j < ne; j++) {
            if (Labels[i] == Labels[j]) {
                throw ParameterException("duplicate Leaf label \"" + Labels[i] + "\", cannot use initial tree");
            }
        }

    if (ne == 3) {
        throw ParameterException("Sorry, this procedure makes no sense for 3 leaves only.");
    }

    /* allocate internal working memory */
    //Edge = new edge_t[2 * ne - 3];
    //BestEdge = new edge_t[2 * ne - 3];
    //EdgeC = new edgec_t[2 * ne - 3];
    Edge = vector<edge_t>(2 * ne - 3);
    BestEdge = vector<edge_t>(2 * ne - 3);
    EdgeC = vector<edgec_t> (2 * ne - 3);
    for (int i = 0; i < 2 * ne - 3; i++)
        EdgeC[i].alloc(NT);

    ShortestLabel = IntVector::Zero(ne);
    inc = IntMatrix::Zero(ne - 2, 3);
    ConShortestPathC = DblMatrix::Zero(2 * ne - 2, NT);
    globminA = DblVector::Zero(NT);
    tmpA = DblVector::Zero(NT);


    if (tree.n_children() != 2) {
        throw ParameterException("Initial tree is not bifurcating");
    } else {
        int si = MapTree2InternalC(tree[0]);
        int sj = MapTree2InternalC(tree[1]);
        if (si >= ne) inc(si - ne, 2) = NewEdgeIndex;
        if (sj >= ne) inc(sj - ne, 2) = NewEdgeIndex;
        EdgeC[NewEdgeIndex].From = si;
        EdgeC[NewEdgeIndex].To = sj;
        double tmp = tree[0].getBranchLength() + tree[1].getBranchLength();
        for (int k = 0; k < NT; k++)
            EdgeC[NewEdgeIndex].len[k] = tmp;
        if (NewInternalNode != ne - 2 || NewEdgeIndex != 2 * ne - 4) {
            throw ParameterException("Initial tree is not of the right size");
        }
    }

    initialized = true;
}

//MinSquareTreeCollection::~MinSquareTreeCollection() {
//    if (EdgeC) {
//        delete[] EdgeC;
//    }
//    if (BestEdge) {
//        delete[] BestEdge;
//    }
//    if (Edge) {
//        delete[] Edge;
//    }
//}

void MinSquareTreeCollection::compute(bool KeepTopology, int iter, bool quiet) throw(RuntimeException) {
    int nswap, i, loop_counter=0;
    double d1, d2, d0;
    IntVector Skip;
    DblMatrix L, M, W;

    if (!initialized) {
        throw RuntimeException("Not initialised.");
    }

    Skip = IntVector::Zero(NT);  // flag for incomplete families
    L = DblMatrix::Zero(NT, 7);  // stores the 7 branch lengths
    M = DblMatrix::Zero(NT, 12); // for 10 distances per family
    W = DblMatrix::Zero(NT, 12); // for 10 weights per family

    IncidencesC();

    globmin = d1 = d0 = DistanceFitCollection();

    if (!quiet) {
        printf("TreeCollection\n");
        printf("--> initial distance fit: %.8g\n", d1);
        printf("--> initial log-likelihood: %.8g\n", LogLikelihoodFitCollection());
    }

    if (KeepTopology) {
        if (!quiet) {
            printf("Beginning branch length optimisation on fixed topology with max %d iterations\n", iter);
        }
        IncidencesC();
        d1 = DBL_MAX;
        i = 4;
        while (i-- > 0) {
            FitLabeledEdgesC(1);
        }
        d2 = DistanceFitCollection();
        i = iter;
        while (i-- > 0 && ReallyLessCollection(d2, d1)) {
            FitLabeledEdgesC(1);
            d1 = d2;
            d2 = DistanceFitCollection();
            if (!quiet) printf("after fitting edges with FitLabeledEdgesC, fit %.12g (d1=%.12g)\n", d2, d1);
        }
        IncidencesC();
        nswap = FiveOptimCollection(1, L, Skip, M, W);
    }
    else {
        i = iter;
        if (!quiet) {
            printf("Beginning topology optimisation with up to %d iterations of branch length optimisation per topology\n", i);
        }
        IncidencesC();

        for (int i=0; i < iter; ++i){
            FitLabeledEdgesC(1);
            d2 = DistanceFitCollection();
            if (ReallyLessCollection(d2, d0))
                d0 = d2;
            else
                break;
        }
        if (!quiet) {
            printf("--> Distance fit   = %.12g\n", d2);
            printf("--> log-likelihood = %.12g\n", LogLikelihoodFitCollection());
        }

        while (i-- > 0) {
            loop_counter += 1;
            if (!quiet) {
                printf("\n\nLOOP %d\n", loop_counter);
                printf("1) Topology optimisation (i = %d)\n", i);
            }
            if (ne > 4) {
                nswap = FiveOptimCollection(0, L, Skip, M, W);
                IncidencesC();
                FitLabeledEdgesC(1);
                d2 = DistanceFitCollection();
                if (!quiet) {
                    printf("--> Distance fit   = %.8g (%d swaps)\n", d2, nswap);
                    printf("--> log-likelihood = %.12g\n", LogLikelihoodFitCollection());
                }
            } else {
                nswap = FourOptimCollection(0, L, Skip, M, W);
                IncidencesC();
                FitLabeledEdgesC(1);
                d2 = DistanceFitCollection();
                if (!quiet) {
                    printf("--> Distance fit   = %.8g (%d swaps)\n", d2, nswap);
                    printf("--> log-likelihood = %.12g\n", LogLikelihoodFitCollection());
                }
            }
            if (!quiet) {
                printf("2) Branch length optimisation\n");
            }
            IncidencesC();
            FitLabeledEdgesC(1);
            d2 = DistanceFitCollection();

            for (int i=0; i < iter; ++i) {
                FitLabeledEdgesC(1);
                d2 = DistanceFitCollection();
                if (!quiet) printf("after fitting edges with FitLabeledEdgesC, %.8g\n", d2);
                if (ReallyLessCollection(d2, d0))
                    d0 = d2;
                else
                    break;
            }

            if (!quiet) {
                printf("--> Distance fit   = %.8g\n", d2);
                printf("--> log-likelihood = %.12g\n", LogLikelihoodFitCollection());
            }
            if (!quiet) {
                printf("3) Check if topology changed\n");
            }
            if (nswap > 0) {
                i = iter;
                if (!quiet) {
                    printf("--> Topology changed: reset i to %d\n", i);
                }
                continue;
            }

            if (!quiet) {
                printf("4) Check for convergence (d1 = %f, d2 = %f)\n   If d2 < d1, continue optimising, else terminate\n", d1,
                       d2);
            }
            if (ReallyLessCollection(d2, d1)) {
                if (!quiet) {
                    printf("--> Continue\n");
                }
                d1 = d2;
                bestTree = getPhyTree();
                bestScore = d2;
                bestLnL = LogLikelihoodFitCollection();
                continue;
            }
            if (!quiet) {
                printf("--> Terminate\n");
            }
            break;
        }
    }
    if (i <= 0)
        printf("Warning: stopped after %d iterations without swap, not yet at minimum\n", iter);

    d1 = ne<=3 ? d1 : d1/(ne-2)/(ne-3)*2;
    MST_Qual = d1;
    computed = true;
}

bool MinSquareTreeCollection::isInitialized() {
    return this->initialized;
}

bool MinSquareTreeCollection::isComputed() {
    return this->computed;
}

double MinSquareTreeCollection::getScore() {
    return bestScore;
}

double MinSquareTreeCollection::getLogLikelihood() {
   return bestLnL;
}

void MinSquareTreeCollection::getTree() {
    PhyTree::TreePtr tree = bestTree;
    if (tree) tree->print(this->newick);
}
