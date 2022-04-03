#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

using Libdl

global podlib

function loadpod(mdppath::String)    
    libexist = 0
    if Sys.isapple()
        libpath = mdppath * "src/Potential/cpuPOD.dylib";
    elseif Sys.isunix()
        libpath = mdppath * "src/Potential/cpuPOD.so";
    elseif Sys.windows()
        libpath = mdppath * "src/Potential/cpuPOD.so";
    end
    if isfile(libpath)
        libexist = 1;
    end        
    if libexist==0            
        cdir = pwd();
        cd(mdppath * "src/Potential")        
        cstr = `g++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuPOD.cpp -o cpuPOD.o`                
        run(cstr)
        if Sys.isapple()            
            cstr = `g++ -shared cpuPOD.o -o cpuPOD.dylib`            
        else
            cstr = `g++ -shared cpuPOD.o -o cpuPOD.so`            
        end 
        run(cstr)
        cd(cdir)                   
    end
    global podlib = Libdl.dlopen(libpath)
end

function closepod()
    Libdl.dlclose(podlib)
end

function podneighborlist(ai, aj, numneigh, y, rcutsq, nx::Int32, N::Int32, dim::Int32)
    count = ccall(Libdl.dlsym(podlib, :neighborlist), Cint, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
            Ptr{Cdouble}, Cdouble, Cint, Cint, Cint), ai, aj, numneigh, y, rcutsq, nx, N, dim)
    return count
end

function podtally3(eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
                    nrbf::Int32, nabf::Int32, natom::Int32, N::Int32)
    ccall(Libdl.dlsym(podlib, :podtally3), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
            eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
            nrbf, nabf, natom, N)    
end

# void podtally3(double *eatom, double *fatom, double *vatom, double *xij, double *xik, double *uij, 
#              double *uik, double *uijk, double *wij, double *wik, double *wijk, int *ai, int *aj,
#              int *ak, int nrbf, int nabf, int natom, int N)