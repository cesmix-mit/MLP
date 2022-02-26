/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

// clang++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuPOD.cpp -o cpuPOD.o
// clang++ -shared cpuPOD.o -o cpuPOD.dylib (MacOS system)
// clang++ -shared cpuPOD.o -o cpuPOD.so (Linux system)

#ifndef CPUPOD
#define CPUPOD

#include "cpuPOD.h"

#include <stdio.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

int neighborlist(int *ai, int *aj, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
{
    int k = 0;
    for (int i = 0; i<nx; i++) {
        double *ri = &r[i*dim];
        int inc = 0;
        for (int j=0; j<N; j++) {
            double *rj = &r[dim*j];                        
            double rijsq = (ri[0]-rj[0])*(ri[0]-rj[0]) + (ri[1]-rj[1])*(ri[1]-rj[1]) + (ri[2]-rj[2])*((ri[2]-rj[2]));
            if  ((rijsq > 1e-12) && (rijsq <= rcutsq))  { 
                inc += 1;                                
                ai[k] = i+1;
                aj[k] = j+1;          
                k += 1;                                                  
            }
        }
        numneigh[i] = inc; 
    }
    return k; 
}

void podtally3(double *eatom, double *fatom, double *vatom, double *xij, double *xik, double *uij, 
             double *uik, double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj,
             int *ak, int nrbf, int nabf, int natom, int N)
{
    int K = -1;
    for (int m =0; m<nrbf; m++) {         
        K += 1;
        double *uj = &uij[N*m];
        double *uk = &uik[N*m];
        double *wj = &wij[3*N*m];
        double *wk = &wik[3*N*m];
        for (int n=0; n<N; n++) {                    
            int k = ai[n] + natom*K;            
            eatom[k] += uj[n]*uk[n];
            double *wj3 = &wj[3*n];
            double *wk3 = &wk[3*n];
            double fij0 = wj3[0]*uk[n];
            double fij1 = wj3[1]*uk[n];
            double fij2 = wj3[2]*uk[n];
            double fik0 = uj[n]*wk3[0];
            double fik1 = uj[n]*wk3[1];
            double fik2 = uj[n]*wk3[2];            
            // fatom[k0] += fij0 + fik0;
            // fatom[k1] += fij1 + fik1;
            // fatom[k2] += fij2 + fik2;
        }
        for (int p=0; p<nabf; p++) {               
            K = K + 1;                
            double *u = &uijk[N*p];
            double *w = &wijk[6*N*p];
            for (int n=0; n<N; n++) {      
                int k = ai[n] + natom*K;                   
                double ujk = uj[n]*uk[n];
                double ujn = uj[n]*u[n];
                double ukn = uk[n]*u[n];
                eatom[k] += ujk*u[n];             
                double *wj3 = &wj[3*n];
                double *wk3 = &wk[3*n];
                double *w6 = &w[6*n];
                double fij0 = wj3[0]*ukn + ujk*w6[0];
                double fij1 = wj3[1]*ukn + ujk*w6[1];
                double fij2 = wj3[2]*ukn + ujk*w6[2];
                double fik0 = wk3[0]*ujn + ujk*w6[3];
                double fik1 = wk3[1]*ujn + ujk*w6[4];
                double fik2 = wk3[2]*ujn + ujk*w6[5];            
            }
        }
    }

    // double fij0, fij1, fij2, xij0, xij1, xij2;
    // double fik0, fik1, fik2, xik0, xik1, xik2;
    // double v0, v1, v2, v3, v4, v5;
    // double eijk, uj, uk, ujk, wij0, wij1, wij2, wik0, wik1, wik2;
            
    // int K = -1;
    // for (int m =0; m<nrbf; m++) {                    
    //     K += 1;
    //     int Nm3 = 3*N*m;        
    //     int nK3 = 3*natom*K;
    //     int nK6 = 6*natom*K;
    //     for (int n=0; n<N; n++) {        
    //         int nm = n + N*m;
    //         // int nm0 = 0 + 3*n + Nm3;
    //         // int nm1 = 1 + 3*n + Nm3;
    //         // int nm2 = 2 + 3*n + Nm3;
    //         // int n0 = 0 + 3*n;
    //         // int n1 = 1 + 3*n;
    //         // int n2 = 2 + 3*n;

    //         uj = uij[nm];
    //         uk = uik[nm];
    //         ujk = uj*uk;
    //         double *xij3 = &xij[3*n];
    //         double *xik3 = &xik[3*n];
    //         double *wij3 = &wij[3*nm];
    //         double *wik3 = &wik[3*nm];
    //         xij0 = xij3[0];
    //         xij1 = xij3[1];
    //         xij2 = xij3[2];
    //         xik0 = xik3[0];
    //         xik1 = xik3[1];
    //         xik2 = xik3[2];
    //         wij0 = wij3[0];
    //         wij1 = wij3[1];
    //         wij2 = wij3[2];
    //         wik0 = wik3[0];
    //         wik1 = wik3[1];
    //         wik2 = wik3[2];
    //         // wij0 = wij[nm0];
    //         // wij1 = wij[nm1];
    //         // wij2 = wij[nm2];
    //         // wik0 = wik[nm0];
    //         // wik1 = wik[nm1];
    //         // wik2 = wik[nm2];

    //         eijk = uj*uk;
    //         fij0 = wij0*uk;
    //         fij1 = wij1*uk;
    //         fij2 = wij2*uk;
    //         fik0 = uj*wik0;
    //         fik1 = uj*wik1;
    //         fik2 = uj*wik2;                        

    //         v0 = xij0*fij0 + xik0*fik0;
    //         v1 = xij1*fij1 + xik1*fik1;        
    //         v2 = xij2*fij2 + xik2*fik2;
    //         v3 = xij1*fij2 + xik1*fik2;
    //         v4 = xij0*fij2 + xik0*fik2;            
    //         v5 = xij0*fij1 + xik0*fik1;                                

    //         int i1 = ai[n];
    //         int k = i1 + natom*K;
    //         int k0 = 0 + 3*i1 + nK3;
    //         int k1 = 1 + 3*i1 + nK3;
    //         int k2 = 2 + 3*i1 + nK3;
    //         int l0 = 0 + 6*i1 + nK6;
    //         int l1 = 1 + 6*i1 + nK6;
    //         int l2 = 2 + 6*i1 + nK6;
    //         int l3 = 3 + 6*i1 + nK6;
    //         int l4 = 4 + 6*i1 + nK6;
    //         int l5 = 5 + 6*i1 + nK6;

    //         eatom[k] += eijk;
    //         // fatom[k0] += fij0 + fik0;
    //         // fatom[k1] += fij1 + fik1;
    //         // fatom[k2] += fij2 + fik2;
    //         // vatom[l0] += v0; 
    //         // vatom[l1] += v1;
    //         // vatom[l2] += v2; 
    //         // vatom[l3] += v3;
    //         // vatom[l4] += v4; 
    //         // vatom[l5] += v5;        
            
    //         // i1 = aj[n];
    //         // k = i1 + natom*K;
    //         // k0 = 0 + 3*i1 + nK3;
    //         // k1 = 1 + 3*i1 + nK3;
    //         // k2 = 2 + 3*i1 + nK3;
    //         // l0 = 0 + 6*i1 + nK6;
    //         // l1 = 1 + 6*i1 + nK6;
    //         // l2 = 2 + 6*i1 + nK6;
    //         // l3 = 3 + 6*i1 + nK6;
    //         // l4 = 4 + 6*i1 + nK6;
    //         // l5 = 5 + 6*i1 + nK6;
    //         // fatom[k0] -= fij0;
    //         // fatom[k1] -= fij1;
    //         // fatom[k2] -= fij2;
    //         // vatom[l0] += v0; 
    //         // vatom[l1] += v1;
    //         // vatom[l2] += v2; 
    //         // vatom[l3] += v3;
    //         // vatom[l4] += v4; 
    //         // vatom[l5] += v5;        

    //         // i1 = ak[n];
    //         // k = i1 + natom*K;
    //         // k0 = 0 + 3*i1 + nK3;
    //         // k1 = 1 + 3*i1 + nK3;
    //         // k2 = 2 + 3*i1 + nK3;
    //         // l0 = 0 + 6*i1 + nK6;
    //         // l1 = 1 + 6*i1 + nK6;
    //         // l2 = 2 + 6*i1 + nK6;
    //         // l3 = 3 + 6*i1 + nK6;
    //         // l4 = 4 + 6*i1 + nK6;
    //         // l5 = 5 + 6*i1 + nK6;
    //         // fatom[k0] -= fik0;   
    //         // fatom[k1] -= fik1;   
    //         // fatom[k2] -= fik2;   
    //         // vatom[l0] += v0; 
    //         // vatom[l1] += v1;
    //         // vatom[l2] += v2; 
    //         // vatom[l3] += v3;
    //         // vatom[l4] += v4; 
    //         // vatom[l5] += v5;        
    //     }

    //     for (int p=0; p<nabf; p++) {               
    //         K = K + 1;                
    //         nK3 = 3*natom*K;
    //         nK6 = 6*natom*K;
    //         for (int n=0; n<N; n++) {       
    //             int nm = n + N*m;
    //             // int nm0 = 0 + 3*n + Nm3;
    //             // int nm1 = 1 + 3*n + Nm3;
    //             // int nm2 = 2 + 3*n + Nm3;
    //             // int n0 = 0 + 3*n;
    //             // int n1 = 1 + 3*n;
    //             // int n2 = 2 + 3*n;
    //             int np = n + N*p;
    //             int np0 = 0 + 6*n + 6*N*p;
    //             int np1 = 1 + 6*n + 6*N*p;
    //             int np2 = 2 + 6*n + 6*N*p;
    //             int np3 = 3 + 6*n + 6*N*p;
    //             int np4 = 4 + 6*n + 6*N*p;
    //             int np5 = 5 + 6*n + 6*N*p;      

    //             uj = uij[nm];
    //             uk = uik[nm];
    //             ujk = uj*uk;
    //             double *xij3 = &xij[3*n];
    //             double *xik3 = &xik[3*n];
    //             double *wij3 = &wij[3*nm];
    //             double *wik3 = &wik[3*nm];
    //             xij0 = xij3[0];
    //             xij1 = xij3[1];
    //             xij2 = xij3[2];
    //             xik0 = xik3[0];
    //             xik1 = xik3[1];
    //             xik2 = xik3[2];
    //             wij0 = wij3[0];
    //             wij1 = wij3[1];
    //             wij2 = wij3[2];
    //             wik0 = wik3[0];
    //             wik1 = wik3[1];
    //             wik2 = wik3[2];
    //             // wij0 = wij[nm0];
    //             // wij1 = wij[nm1];
    //             // wij2 = wij[nm2];
    //             // wik0 = wik[nm0];
    //             // wik1 = wik[nm1];
    //             // wik2 = wik[nm2];

    //             double u = uijk[np];   
    //             eijk = ujk*u;                
    //             fij0 = wij0*uk*u + ujk*wijk[np0];
    //             fij1 = wij1*uk*u + ujk*wijk[np1];
    //             fij2 = wij2*uk*u + ujk*wijk[np2];
    //             fik0 = uj*wik0*u + ujk*wijk[np3];
    //             fik1 = uj*wik1*u + ujk*wijk[np4];
    //             fik2 = uj*wik2*u + ujk*wijk[np5];
           
    //             v0 = xij0*fij0 + xik0*fik0;
    //             v1 = xij1*fij1 + xik1*fik1;        
    //             v2 = xij2*fij2 + xik2*fik2;
    //             v3 = xij1*fij2 + xik1*fik2;
    //             v4 = xij0*fij2 + xik0*fik2;            
    //             v5 = xij0*fij1 + xik0*fik1;                                

    //             int i1 = ai[n];
    //             int k = i1 + natom*K;
    //             int k0 = 0 + 3*i1 + nK3;
    //             int k1 = 1 + 3*i1 + nK3;
    //             int k2 = 2 + 3*i1 + nK3;
    //             int l0 = 0 + 6*i1 + nK6;
    //             int l1 = 1 + 6*i1 + nK6;
    //             int l2 = 2 + 6*i1 + nK6;
    //             int l3 = 3 + 6*i1 + nK6;
    //             int l4 = 4 + 6*i1 + nK6;
    //             int l5 = 5 + 6*i1 + nK6;
    //             eatom[k] += eijk;
    //             // fatom[k0] += fij0 + fik0;
    //             // fatom[k1] += fij1 + fik1;
    //             // fatom[k2] += fij2 + fik2;
    //             // vatom[l0] += v0; 
    //             // vatom[l1] += v1;
    //             // vatom[l2] += v2; 
    //             // vatom[l3] += v3;
    //             // vatom[l4] += v4; 
    //             // vatom[l5] += v5;        
                
    //             // i1 = aj[n];
    //             // k = i1 + natom*K;
    //             // k0 = 0 + 3*i1 + nK3;
    //             // k1 = 1 + 3*i1 + nK3;
    //             // k2 = 2 + 3*i1 + nK3;
    //             // l0 = 0 + 6*i1 + nK6;
    //             // l1 = 1 + 6*i1 + nK6;
    //             // l2 = 2 + 6*i1 + nK6;
    //             // l3 = 3 + 6*i1 + nK6;
    //             // l4 = 4 + 6*i1 + nK6;
    //             // l5 = 5 + 6*i1 + nK6;
    //             // fatom[k0] -= fij0;
    //             // fatom[k1] -= fij1;
    //             // fatom[k2] -= fij2;
    //             // vatom[l0] += v0; 
    //             // vatom[l1] += v1;
    //             // vatom[l2] += v2; 
    //             // vatom[l3] += v3;
    //             // vatom[l4] += v4; 
    //             // vatom[l5] += v5;      

    //             // i1 = ak[n];
    //             // k = i1 + natom*K;
    //             // k0 = 0 + 3*i1 + nK3;
    //             // k1 = 1 + 3*i1 + nK3;
    //             // k2 = 2 + 3*i1 + nK3;
    //             // l0 = 0 + 6*i1 + nK6;
    //             // l1 = 1 + 6*i1 + nK6;
    //             // l2 = 2 + 6*i1 + nK6;
    //             // l3 = 3 + 6*i1 + nK6;
    //             // l4 = 4 + 6*i1 + nK6;
    //             // l5 = 5 + 6*i1 + nK6;
    //             // fatom[k0] -= fik0;   
    //             // fatom[k1] -= fik1;   
    //             // fatom[k2] -= fik2;   
    //             // vatom[l0] += v0; 
    //             // vatom[l1] += v1;
    //             // vatom[l2] += v2; 
    //             // vatom[l3] += v3;
    //             // vatom[l4] += v4; 
    //             // vatom[l5] += v5;            
    //        }
    //     }
    // }
    // for (int n=0; n<N; n++) {
    //     int K = -1;
    //     for (int m =0; m<nrbf; m++) {
    //         K += 1;
    //         int nm = n + N*m;
    //         int nm0 = 0 + 3*n + 3*N*m;
    //         int nm1 = 1 + 3*n + 3*N*m;
    //         int nm2 = 2 + 3*n + 3*N*m;
    //         int n0 = 0 + 3*n;
    //         int n1 = 1 + 3*n;
    //         int n2 = 2 + 3*n;

    //         double uj = uij[nm];
    //         double uk = uik[nm];
    //         double ujk = uj*uk;
    //         double wij0 = wij[nm0];
    //         double wij1 = wij[nm1];
    //         double wij2 = wij[nm2];
    //         double wik0 = wik[nm0];
    //         double wik1 = wik[nm1];
    //         double wik2 = wik[nm2];

    //         double eijk = uj*uk;
    //         fij[0] = wij0*uk;
    //         fij[1] = wij1*uk;
    //         fij[2] = wij2*uk;
    //         fik[0] = uj*wik0;
    //         fik[1] = uj*wik1;
    //         fik[2] = uj*wik2;            

    //         double v0 = xij[n0]*fij[0] + xik[n0]*fik[0];
    //         double v1 = xij[n1]*fij[1] + xik[n1]*fik[1];        
    //         double v2 = xij[n2]*fij[2] + xik[n2]*fik[2];
    //         double v3 = xij[n1]*fij[2] + xik[n1]*fik[2];
    //         double v4 = xij[n0]*fij[2] + xik[n0]*fik[2];            
    //         double v5 = xij[n0]*fij[1] + xik[n0]*fik[1];                                

    //         int i1 = ai[n];
    //         int k = i1 + natom*K;
    //         int k0 = 0 + 3*i1 + 3*natom*K;
    //         int k1 = 1 + 3*i1 + 3*natom*K;
    //         int k2 = 2 + 3*i1 + 3*natom*K;
    //         int l0 = 0 + 6*i1 + 6*natom*K;
    //         int l1 = 1 + 6*i1 + 6*natom*K;
    //         int l2 = 2 + 6*i1 + 6*natom*K;
    //         int l3 = 3 + 6*i1 + 6*natom*K;
    //         int l4 = 4 + 6*i1 + 6*natom*K;
    //         int l5 = 5 + 6*i1 + 6*natom*K;

    //         eatom[k] += eijk;
    //         fatom[k0] += fij[0] + fik[0];
    //         fatom[k1] += fij[1] + fik[1];
    //         fatom[k2] += fij[2] + fik[2];
    //         vatom[l0] += v0; 
    //         vatom[l1] += v1;
    //         vatom[l2] += v2; 
    //         vatom[l3] += v3;
    //         vatom[l4] += v4; 
    //         vatom[l5] += v5;        
            
    //         i1 = aj[n];
    //         k = i1 + natom*K;
    //         k0 = 0 + 3*i1 + 3*natom*K;
    //         k1 = 1 + 3*i1 + 3*natom*K;
    //         k2 = 2 + 3*i1 + 3*natom*K;
    //         l0 = 0 + 6*i1 + 6*natom*K;
    //         l1 = 1 + 6*i1 + 6*natom*K;
    //         l2 = 2 + 6*i1 + 6*natom*K;
    //         l3 = 3 + 6*i1 + 6*natom*K;
    //         l4 = 4 + 6*i1 + 6*natom*K;
    //         l5 = 5 + 6*i1 + 6*natom*K;
    //         fatom[k0] -= fij[0];
    //         fatom[k1] -= fij[1];
    //         fatom[k2] -= fij[2];
    //         vatom[l0] += v0; 
    //         vatom[l1] += v1;
    //         vatom[l2] += v2; 
    //         vatom[l3] += v3;
    //         vatom[l4] += v4; 
    //         vatom[l5] += v5;        

    //         i1 = ak[n];
    //         k = i1 + natom*K;
    //         k0 = 0 + 3*i1 + 3*natom*K;
    //         k1 = 1 + 3*i1 + 3*natom*K;
    //         k2 = 2 + 3*i1 + 3*natom*K;
    //         l0 = 0 + 6*i1 + 6*natom*K;
    //         l1 = 1 + 6*i1 + 6*natom*K;
    //         l2 = 2 + 6*i1 + 6*natom*K;
    //         l3 = 3 + 6*i1 + 6*natom*K;
    //         l4 = 4 + 6*i1 + 6*natom*K;
    //         l5 = 5 + 6*i1 + 6*natom*K;
    //         fatom[k0] -= fik[0];   
    //         fatom[k1] -= fik[1];   
    //         fatom[k2] -= fik[2];   
    //         vatom[l0] += v0; 
    //         vatom[l1] += v1;
    //         vatom[l2] += v2; 
    //         vatom[l3] += v3;
    //         vatom[l4] += v4; 
    //         vatom[l5] += v5;        

    //         for (int p=0; p<nabf; p++) {               
    //             K = K + 1;                
    //             int np = n + N*p;
    //             int np0 = 0 + 6*n + 6*N*p;
    //             int np1 = 1 + 6*n + 6*N*p;
    //             int np2 = 2 + 6*n + 6*N*p;
    //             int np3 = 3 + 6*n + 6*N*p;
    //             int np4 = 4 + 6*n + 6*N*p;
    //             int np5 = 5 + 6*n + 6*N*p;            
    //             double u = uijk[np];   
    //             eijk = ujk*u;                
    //             fij[0] = wij0*uk*u + ujk*wijk[np0];
    //             fij[1] = wij1*uk*u + ujk*wijk[np1];
    //             fij[2] = wij2*uk*u + ujk*wijk[np2];
    //             fik[0] = uj*wik0*u + ujk*wijk[np3];
    //             fik[1] = uj*wik1*u + ujk*wijk[np4];
    //             fik[2] = uj*wik2*u + ujk*wijk[np5];
           
    //             v0 = xij[n0]*fij[0] + xik[n0]*fik[0];
    //             v1 = xij[n1]*fij[1] + xik[n1]*fik[1];        
    //             v2 = xij[n2]*fij[2] + xik[n2]*fik[2];
    //             v3 = xij[n1]*fij[2] + xik[n1]*fik[2];
    //             v4 = xij[n0]*fij[2] + xik[n0]*fik[2];            
    //             v5 = xij[n0]*fij[1] + xik[n0]*fik[1];                                

    //             i1 = ai[n];
    //             k = i1 + natom*K;
    //             k0 = 0 + 3*i1 + 3*natom*K;
    //             k1 = 1 + 3*i1 + 3*natom*K;
    //             k2 = 2 + 3*i1 + 3*natom*K;
    //             l0 = 0 + 6*i1 + 6*natom*K;
    //             l1 = 1 + 6*i1 + 6*natom*K;
    //             l2 = 2 + 6*i1 + 6*natom*K;
    //             l3 = 3 + 6*i1 + 6*natom*K;
    //             l4 = 4 + 6*i1 + 6*natom*K;
    //             l5 = 5 + 6*i1 + 6*natom*K;
    //             eatom[k] += eijk;
    //             fatom[k0] += fij[0] + fik[0];
    //             fatom[k1] += fij[1] + fik[1];
    //             fatom[k2] += fij[2] + fik[2];
    //             vatom[l0] += v0; 
    //             vatom[l1] += v1;
    //             vatom[l2] += v2; 
    //             vatom[l3] += v3;
    //             vatom[l4] += v4; 
    //             vatom[l5] += v5;        
                
    //             i1 = aj[n];
    //             k = i1 + natom*K;
    //             k0 = 0 + 3*i1 + 3*natom*K;
    //             k1 = 1 + 3*i1 + 3*natom*K;
    //             k2 = 2 + 3*i1 + 3*natom*K;
    //             l0 = 0 + 6*i1 + 6*natom*K;
    //             l1 = 1 + 6*i1 + 6*natom*K;
    //             l2 = 2 + 6*i1 + 6*natom*K;
    //             l3 = 3 + 6*i1 + 6*natom*K;
    //             l4 = 4 + 6*i1 + 6*natom*K;
    //             l5 = 5 + 6*i1 + 6*natom*K;
    //             fatom[k0] -= fij[0];
    //             fatom[k1] -= fij[1];
    //             fatom[k2] -= fij[2];
    //             vatom[l0] += v0; 
    //             vatom[l1] += v1;
    //             vatom[l2] += v2; 
    //             vatom[l3] += v3;
    //             vatom[l4] += v4; 
    //             vatom[l5] += v5;      

    //             i1 = ak[n];
    //             k = i1 + natom*K;
    //             k0 = 0 + 3*i1 + 3*natom*K;
    //             k1 = 1 + 3*i1 + 3*natom*K;
    //             k2 = 2 + 3*i1 + 3*natom*K;
    //             l0 = 0 + 6*i1 + 6*natom*K;
    //             l1 = 1 + 6*i1 + 6*natom*K;
    //             l2 = 2 + 6*i1 + 6*natom*K;
    //             l3 = 3 + 6*i1 + 6*natom*K;
    //             l4 = 4 + 6*i1 + 6*natom*K;
    //             l5 = 5 + 6*i1 + 6*natom*K;
    //             fatom[k0] -= fik[0];   
    //             fatom[k1] -= fik[1];   
    //             fatom[k2] -= fik[2];   
    //             vatom[l0] += v0; 
    //             vatom[l1] += v1;
    //             vatom[l2] += v2; 
    //             vatom[l3] += v3;
    //             vatom[l4] += v4; 
    //             vatom[l5] += v5;            
    //        }
    //     }
    // }
}

#endif

