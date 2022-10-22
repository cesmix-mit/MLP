#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int neighborlist(int *ai, int *aj, int *numneigh, double *r, double rcutsq, int nx, int N, int dim);

int podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, 
        int *numneighsum, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);

void pod3desc(double *eatom, double *fatom, double *y, double *x, double *a1, double *a2, double *a3,
             double *phi, double *gamma0, double *tmpmem, double rin, double rcut, int *tmpint, 
            int *elemind, int *pairlist, int *pairnum, int *pairnumsum, int *atomtype, int *alist, 
            int *pbc, int *pdegree, int ngm, int nrbf, int nabf, int nelements, int natom);

void pod3bodydesc(double *eatom, double *fatom, double *x, double *a1, double *a2, double *a3,
             double *phi, double *gamma0, double rin, double rcut, int *atomtype, 
            int *pbc, int *pdegree, int ngm, int nrbf, int nabf, int nelements, int natom);

void poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
            double *fatom3, double *x, double *a1, double *a2, double *a3, double *phi2, double *phi3, 
            double *gamma0, double rin, double rcut, int *atomtype, int *pbc, int *pdegree2, 
            int *pdegree3, int ngm, int nrbf2, int nrbf3, int nabf, int nelements, int natom);

void bispectrumderiv(double *dbsa, double *Ar, double *Ai, double *dAr, double *dAi, double *cg, 
    int *pairnum, int *indl, int *indm, int *rowm, int K, int J, int Q, int M, int N, int Nij);

int cpuNeighPairList(int *pairnum, int *pairnumsum, int *pairlist, double *x, double rcutsq, 
        int *neighlist,  int *neighnumsum, int inum,  int dim);

void cpuNeighPairs(double *xij, double *x, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *atomtype, int *alist, int inum, int dim);

void cpuNeighTripletList(int *tripletlist, int *tripletnum, int *tripletnumsum, int *pairlist, 
     int *pairnum, int *pairnumsum, int inum);

void cpuNeighTriplets(double *xij, double *xik, double *x, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletlist, int *tripletnumsum, 
      int *alist,  int *atomtype, int inum, int dim);

void podhybrid23(double* d23, double *dd23, double* d2, double *d3, double* dd2, double *dd3, 
        int M2, int M3, int N);
        
void makeindjk(int *indj, int *indk, int n);

void makejk(double *uij, double *uik, double *wij, double *wik, double *e2ij, double *f2ij, 
        double *ei, double *f1, double *f2, double *f3, int *pairnum, int *tripletnum, 
        int *indj, int *indk, int Nj, int M, int inum, int Nij, int Nijk);

void radialbasis0(double *rbf, double *xij, double *scalefac,
                double rin, double rmax, int pdegree, int K, int M, int N);

void radialbasis(double *rbf, double *drbf, double *xij, double *scalefac,
                double rin, double rmax, int pdegree, int K, int M, int N);

void cosinbasis(double *abf, double *dabf, double *xij, double *xik, int nabf, int N);

void podtally2(double *eatom, double *fatom, double *vatom, double *rij, double *eij, double *fij, 
               int *ai, int *aj, int *ti, int *tj, int *ind, int S, int natom, int M, int N);

void podtally2b(double *eatom, double *fatom, double *eij, double *fij, 
               int *ai, int *aj, int *ti, int *tj, int *ind, int S, int natom, int M, int N);

void podtally3(double *eatom, double *fatom, double *vatom, double *xij, double *xik, double *uij, 
             double *uik, double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj,
             int *ak, int nrbf, int nabf, int natom, int N);

void podtally3b(double *eatom, double *fatom, double *uij, double *uik, 
                double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj, int *ak, 
                int *ti, int *tj, int *tk, int *ind, int nrbf, int nabf, int natom, int N, int S);

void pod3body(double *eatom, double *fatom, double *y, double *phi, double *gamma0, double *tmpmem, 
             double rin, double rcut, int *tmpint, int *ind, int *pairlist, int *pairnum, 
             int *pairnumsum, int *atomtype, int *alist, int *pdegree, int ngm, int nbf, 
             int nrbf, int nabf, int nelements, int natom, int Nj, int Nij, int Nijk);

void pod3body0(double *eatom, double *y, double *phi, double *gamma0, double *tmpmem, 
             double rin, double rcut, int *tmpint, int *ind, int *pairlist, int *pairnum, 
             int *pairnumsum, int *atomtype, int *alist, int *pdegree, int ngm, int nbf, 
             int nrbf, int nabf, int nelements, int natom, int Nj, int Nij, int Nijk);

#ifdef __cplusplus
}
#endif


