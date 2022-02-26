#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int neighborlist(int *ai, int *aj, int *numneigh, double *r, double rcutsq, int nx, int N, int dim);

void podtally3(double *eatom, double *fatom, double *vatom, double *xij, double *xik, double *uij, 
             double *uik, double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj,
             int *ak, int nrbf, int nabf, int natom, int N);

#ifdef __cplusplus
}
#endif


