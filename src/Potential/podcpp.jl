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
    elseif Sys.iswindows()
        libpath = mdppath * "src/Potential/cpuPOD.dll";
    end
    if isfile(libpath)
        libexist = 1;
    end        
    if libexist==0            
        cdir = pwd();
        cd(mdppath * "src/Potential")        
        cstr = `g++ -std=c++11 -O3 -Wall -mtune=intel -ffast-math -Wextra -pedantic -c -fPIC cpuPOD.cpp -o cpuPOD.o`                
        run(cstr)
        if Sys.isapple()            
            cstr = `g++ -shared cpuPOD.o -o cpuPOD.dylib`            
        elseif Sys.isunix()
            cstr = `g++ -shared cpuPOD.o -o cpuPOD.so`      
        elseif Sys.iswindows()
            cstr = `g++ -shared cpuPOD.o -o cpuPOD.dll`    
        end 
        run(cstr)
        cd(cdir)                   
    end
    global podlib = Libdl.dlopen(libpath)
end

function closepod()
    Libdl.dlclose(podlib)
end

function podneighborlist(ai::Vector{Int32}, aj::Vector{Int32}, numneigh::Vector{Int32}, y, rcutsq, nx::Int32, N::Int32, dim::Int32)
    count = ccall(Libdl.dlsym(podlib, :neighborlist), Cint, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
            Ptr{Cdouble}, Cdouble, Cint, Cint, Cint), ai, aj, numneigh, y, rcutsq, nx, N, dim)
    return count
end

function NeighPairList(pairnum::Vector{Int32}, pairnumsum::Vector{Int32}, pairlist::Vector{Int32}, 
    x, rcutsq, neighlist::Vector{Int32}, neighnumsum::Vector{Int32}, inum::Int32, dim::Int32)
    count = ccall(Libdl.dlsym(podlib, :cpuNeighPairList), Cint, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
            Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
            pairnum, pairnumsum, pairlist, x, rcutsq, neighlist, neighnumsum, inum, dim)
    return count
end

function NeighPairs(xij, x, ai::Vector{Int32}, aj::Vector{Int32}, ti::Vector{Int32}, tj::Vector{Int32}, 
    pairlist::Vector{Int32}, pairnumsum::Vector{Int32}, atomtype::Vector{Int32}, alist::Vector{Int32}, 
    inum::Int32, dim::Int32)
    ccall(Libdl.dlsym(podlib, :cpuNeighPairs), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, 
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
            xij, x, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, alist, inum, dim)    
end

function NeighTripletList(tripletlist::Vector{Int32}, tripletnum::Vector{Int32}, tripletnumsum::Vector{Int32}, 
    pairlist::Vector{Int32}, pairnum::Vector{Int32}, pairnumsum::Vector{Int32},  inum::Int32)
    ccall(Libdl.dlsym(podlib, :cpuNeighTripletList), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint), tripletlist, tripletnum, 
            tripletnumsum, pairlist, pairnum, pairnumsum, inum)    
end

function NeighTriplets(xij, xik, x, ai::Vector{Int32}, aj::Vector{Int32}, ak::Vector{Int32}, 
    ti::Vector{Int32}, tj::Vector{Int32}, tk::Vector{Int32}, tripletlist::Vector{Int32}, 
    tripletnumsum::Vector{Int32}, alist::Vector{Int32}, atomtype::Vector{Int32}, inum::Int32, dim::Int32)
    ccall(Libdl.dlsym(podlib, :cpuNeighTriplets), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint), xij, xik, x, ai, aj, ak, ti, tj, tk, 
            tripletlist, tripletnumsum, alist, atomtype, inum, dim)    
end

# void bispectrumderiv(double *dbsa, double *Ar, double *Ai, double *dAr, double *dAi, double *cg, 
#     int *pairnum, int *indl, int *indm, int *rowm, int K, int J, int Q, int M, int N, int Nij);
function bispectrumderiv(dbsa, Ar, Ai, dAr, dAi, cg, pairnum::Vector{Int32}, indl::Vector{Int32}, 
    indm::Vector{Int32}, rowm::Vector{Int32}, K::Int32, J::Int32, Q::Int32, M::Int32, N::Int32, Nij::Int32)
    ccall(Libdl.dlsym(podlib, :bispectrumderiv), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Cint, Cint, Cint, Cint, Cint, Cint), dbsa, Ar, Ai, dAr, dAi, cg, pairnum, indl, indm, rowm, 
        K, J, Q, M, N, Nij)    
end

function makeindjk(indj::Vector{Int32}, indk::Vector{Int32}, n::Int32)
    ccall(Libdl.dlsym(podlib, :makeindjk), Cvoid, (Ptr{Cint}, Ptr{Cint}, Cint), indj, indk, n)    
end

function makejk(uij, uik, wij, wik, e2ij, f2ej, ei, f1, f2, f3, pairnum::Vector{Int32},
    tripletnum::Vector{Int32}, indj::Vector{Int32}, indk::Vector{Int32}, Nj::Int32, M::Int32, inum::Int32, Nij::Int32, Nijk::Int32)
    ccall(Libdl.dlsym(podlib, :makejk), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},  Cint, Cint, Cint, Cint, Cint), 
    uij, uik, wij, wik, e2ij, f2ej, ei, f1, f2, f3, pairnum, tripletnum, indj, indk, 
    Nj, M, inum, Nij, Nijk)    
end

function radialbasis0(rbf, xij, scalefac, rin, rmax, pdegree::Int32, K::Int32, M::Int32, N::Int32)
    ccall(Libdl.dlsym(podlib, :radialbasis0), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Cdouble, Cdouble, Cint, Cint, Cint, Cint), rbf, xij, scalefac, 
        rin, rmax, pdegree, K, M, N)    
end

function radialbasis(rbf, drbf, xij, scalefac, rin, rmax, pdegree::Int32, K::Int32, M::Int32, N::Int32)
    ccall(Libdl.dlsym(podlib, :radialbasis), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Cdouble, Cdouble, Cint, Cint, Cint, Cint), rbf, drbf, xij, scalefac, 
        rin, rmax, pdegree, K, M, N)    
end

function cosinbasis(abf, dabf, xij, xik, nabf::Int32, N::Int32)
    ccall(Libdl.dlsym(podlib, :cosinbasis), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Cint, Cint), abf, dabf, xij, xik, nabf, N)    
end

function podhybrid23(d23, dd23, d2, d3, dd2, dd3, M2::Int32, M3::Int32, N::Int32)
    ccall(Libdl.dlsym(podlib, :podhybrid23), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint), d23, dd23, d2, d3, dd2, dd3, M2, M3, N)    
end

function podtally3(eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai::Vector{Int32}, 
    aj::Vector{Int32}, ak::Vector{Int32}, nrbf::Int32, nabf::Int32, natom::Int32, N::Int32)
    ccall(Libdl.dlsym(podlib, :podtally3), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
            eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
            nrbf, nabf, natom, N)    
end

function podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai::Vector{Int32}, 
    aj::Vector{Int32}, ak::Vector{Int32}, ti::Vector{Int32}, tj::Vector{Int32}, tk::Vector{Int32}, 
    ind::Vector{Int32}, nrbf::Int32, nabf::Int32, natom::Int32, N::Int32, S::Int32)
    ccall(Libdl.dlsym(podlib, :podtally3b), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
    Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint), eatom, fatom, xij, xik, uij, uik, 
    uijk, wij, wik, wijk, ai, aj, ak, ti, tj, tk, ind, nrbf, nabf, natom, N, S)    
end

function pod3body(eatom, fatom, y, phi, gamma0, tmpmem, rin, rcut, tmpint::Vector{Int32}, 
    ind::Vector{Int32}, pairlist::Vector{Int32}, pairnum::Vector{Int32}, pairnumsum::Vector{Int32}, 
    atomtype::Vector{Int32}, alist::Vector{Int32}, pdegree::Vector{Int32}, ngm::Int32, nbf::Int32, 
    nrbf::Int32, nabf::Int32, nelements::Int32, natom::Int32, Nj::Int32, Nij::Int32, Nijk::Int32)
    ccall(Libdl.dlsym(podlib, :pod3body), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble,
    Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
    Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), eatom, fatom, y, phi, gamma0, 
    tmpmem, rin, rcut, tmpint, ind, pairlist, pairnum, pairnumsum, atomtype, alist, pdegree, 
    ngm, nbf, nrbf, nabf, nelements, natom, Nj, Nij, Nijk)    
end

function pod3body0(eatom, y, phi, gamma0, tmpmem, rin, rcut, tmpint::Vector{Int32}, 
    ind::Vector{Int32}, pairlist::Vector{Int32}, pairnum::Vector{Int32}, pairnumsum::Vector{Int32}, 
    atomtype::Vector{Int32}, alist::Vector{Int32}, pdegree::Vector{Int32}, ngm::Int32, nbf::Int32, 
    nrbf::Int32, nabf::Int32, nelements::Int32, natom::Int32, Nj::Int32, Nij::Int32, Nijk::Int32)
    ccall(Libdl.dlsym(podlib, :pod3body0), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
    Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
    eatom, y, phi, gamma0, tmpmem, rin, rcut, tmpint, ind, pairlist, pairnum, pairnumsum, atomtype, 
    alist, pdegree, ngm, nbf, nrbf, nabf, nelements, natom, Nj, Nij, Nijk)    
end
