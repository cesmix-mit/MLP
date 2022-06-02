/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

// clang++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuPOD.cpp -o cpuPOD.o
// clang++ -shared cpuPOD.o -o cpuPOD.dylib (MacOS system)
// clang++ -shared cpuPOD.o -o cpuPOD.so (Linux system)

#ifndef CPUACG
#define CPUACG

#include "cpuACG.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

void similaritymeasure(double* S, double *M,  int l, int n)
{
    for (int i = 0; i<n; i++)
        for (int j = i; j<n; j++)
        {
            double ta = 0.0;
            for (int k=0; k<l; k++)
                ta += M[k + l*i]*M[k + l*j];
                
        }
}

    // for i = 1:n
    //     display(i)
    //     for j = i:n         
    //         ta = 0.0;
    //         #tb = 0.0;
    //         for k = 1:l                
    //             ta += M[k,i]*M[k,j]; 
    //             #tb += (M[k,i]-M[k,j])*(M[k,i]-M[k,j]);
    //         end
    //         tm = abs(ta);#*exp(-tb);            
    //         S[i,j] = tm;
    //         S[j,i] = tm; 
    //     end
    // end

#endif

