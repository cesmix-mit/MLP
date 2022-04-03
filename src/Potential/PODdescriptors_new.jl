include("sphericalharmonics.jl")
include("clebschgordan.jl")
#include("bispectrum.jl")

@kwdef mutable struct PODdesc 
    name::String = "POD"
    nbody::Int64 = 2
    species::Vector{Symbol} = [:X]
    atomtypes::Vector{Int64} = [1, 1, 1, 1] # a list of atom types 
    pdegree::Vector{Int64} = [5,12,5] # radial and angular polynomial degrees 
    nbasis::Vector{Int64} = [5,5] # number of radial and angular descriptors 
    rcut::Float64 = 5.0    # cut-off radius 
    rin::Float64 = 1.0     # inner cut-off radius 
    gamma0::Vector{Float64} = [-2, 0.0, 2] # radial parameters
    theta0::Vector{Float64} = [0]    # angular parameters
    Phi::Matrix{Float64} = zeros(1,1) # radial projection matrix 
    Psi::Matrix{Float64} = zeros(1,1) # angular projection matrix 
    Lambda::Vector{Float64} = [0.0]   # radial eigenvalues
    Upsilon::Vector{Float64} = [0.0]  # angular eigenvalues
    sh::SHBasis = SHBasis(1)
    indelem::Vector{Int32} = [1]
    tmpmem::Vector{Float64} = [0.0]
    tmpint::Vector{Int32} = [0]
end

function initPOD(pod::PODdesc)
        
    tol = 1e-6
    if pod.nbody > 1
        pod.Phi, pod.Lambda = radialbasisprojection(pod.pdegree[1], pod.rin, pod.rcut, pod.gamma0, pod.pdegree[2], tol)    
        pod.Phi = pod.Phi[:,1:pod.nbasis[1]]    
    end
    if pod.nbody > 2
        #pod.Psi, pod.Upsilon = fourierbasisprojection(pod.pdegree[3], tol)
        pod.Psi, pod.Upsilon = cosinbasisprojection(pod.pdegree[3], pod.theta0, tol)
        pod.Psi = pod.Psi[:,1:pod.nbasis[2]]       
    end

    return pod 
end

function POD(; name::String = "POD",
    nbody::Int64 = 2,
    species::Vector{Symbol} = [:X],
    atomtypes::Vector{Int64} = [1, 1, 1, 1], # a list of atom types 
    pdegree::Vector{Int64} = [5,12,5], # radial and angular polynomial degrees 
    nbasis::Vector{Int64} = [5,5], # number of radial and angular descriptors 
    rcut::Float64 = 5.0,    # cut-off radius 
    rin::Float64 = 1.0,     # inner cut-off radius 
    gamma0::Vector{Float64} = [0.0, 2, 4.0], # radial parameters
    theta0::Vector{Float64} = [0.0],    # angular parameters
    Phi::Matrix{Float64} = zeros(1,1), # radial projection matrix 
    Psi::Matrix{Float64} = zeros(1,1), # angular projection matrix 
    Lambda::Vector{Float64} = [0.0],   # radial eigenvalues
    Upsilon::Vector{Float64} = [0.0],  # angular eigenvalues
    )
        
    S = length(species)
    ind = zeros(Int32,S,S)
    k = 0
    for i = 1:S
        for j = i:S
            k = k + 1
            ind[i,j] = k
            ind[j,i] = k
        end
    end
    indelem = Int32.(ind[:])

    if length(nbasis) > 1
        sh = SHBasis(nbasis[2]) 
    else
        sh = SHBasis(1) 
    end

    pod = initPOD(PODdesc(name, nbody, species, atomtypes, pdegree, nbasis, rcut, rin, 
        gamma0, theta0, Phi, Psi, Lambda, Upsilon, sh, indelem, [1.0], [Int32(0)]))

    return pod 
end

function cutoff(dij, rin, rcut, epsil) 

    y = (dij .- rin)*(1.0/(rcut-rin));
    #fcut = exp.(-1.0./sqrt.((1.0 .- y.^2).^2 .+ epsil))/exp(-1.0);
    fcut = exp.(-1.0./sqrt.((1.0 .- y.^3).^2 .+ epsil))/exp(-1.0);
    #fcut = exp.(-1.0./sqrt.((1.0 .- y.^4).^2 .+ epsil))/exp(-1.0);
    
    #fcut = 0.5*(1.0 .+ cos.(pi*y))
    return fcut 
end

function cutoffderiv(dij, rin, rcut, epsil) 

    rmax = rcut - rin
    y = (dij .- rin)/(rcut-rin);

    # fcut = exp.(-1.0./sqrt.((1.0 .- y.^2).^2 .+ epsil))/exp(-1.0);
    # dfcut = ((2.0/(rmax*exp(-1.0)))*y.*exp.(-1.0./(epsil .+ (y.^2 .- 1.0).^2).^(1/2)).*(y.^2 .- 1.0))./(epsil .+ (y.^2 .- 1.0).^2).^(3/2);

    fcut = exp.(-1.0./sqrt.((1.0 .- y.^3).^2 .+ epsil))/exp(-1.0);
    dfcut = ((3.0/(rmax*exp(-1.0)))*(y.^2).*exp.(-1.0./(epsil .+ (y.^3 .- 1.0).^2).^(1/2)).*(y.^3 .- 1.0))./(epsil .+ (y.^3 .- 1.0).^2).^(3/2);

    # fcut = exp.(-1.0./sqrt.((1.0 .- y.^4).^2 .+ epsil))/exp(-1.0);
    # dfcut = ((4.0/(rmax*exp(-1.0)))*(y.^3).*exp.(-1.0./(epsil .+ (y.^4 .- 1.0).^2).^(1/2)).*(y.^4 .- 1.0))./(epsil .+ (y.^4 .- 1.0).^2).^(3/2);
 
    # fcut = 0.5*(1.0 .+ cos.(pi*y))
    # dfcut = -(pi*sin.(pi*y))/(2.0*rmax)
    return fcut, dfcut  
end


function radialprojection(pdegree, dij, rin, rcut, scalefac, K)

    epsil = 1e-6;        
    N = length(dij);
    M = length(scalefac);
    nbf = pdegree*M + K
    rbf = zeros(N, nbf);

    r = dij .- rin
    rmax = rcut - rin
    fcut = cutoff(dij, rin, rcut, epsil) 
    for j = 1:M
        alpha = scalefac[j];    
        if alpha == 0
            alpha = 1e-3;
        end
        x =  (1.0 .- exp.(-alpha*r/rmax))/(1.0-exp(-alpha));

        for i = 1:pdegree    
            a = i*pi;
            rbf[:,i+(j-1)*pdegree] = (sqrt(2.0/(rmax))/i)*fcut.*sin.(a*x)./r;
        end
    end

    # add radial basis
    for i = 1:K
        n = pdegree*M+i;
        rbf[:,n] = fcut./(dij.^i);
    end

    rbfnew, U, s = eigendecomposition(rbf, dij)
    
    return rbf, rbfnew, U, s
end

function radialbasisprojection(pdegree, rin, rcut, scalefac, K, tol)

    dij = LinRange(rin+1e-6, rcut, 2000)
    dij = collect(dij)

    rbf, rbfnew, U, s = radialprojection(pdegree, dij, rin, rcut, scalefac, K)
    # k = numbasis(s, tol)
    # U = U[:,1:k]
    
    return U, s
end

function fourierprojection(pdegree, theta)
    
    N = length(theta);
    nbf = 2*pdegree;
    abf = zeros(N, nbf);        
    for i = 1:pdegree            
        abf[:,i] = cos.(i*theta) 
        abf[:,i+pdegree] = sin.(i*theta) 
    end
    abfnew, U, s = eigendecomposition(abf, theta)

    return abf, abfnew, U, s
end

function fourierbasisprojection(pdegree, tol)

    theta = collect(LinRange(0.0, pi, 2000))     
    rbf, rbfnew, U, s = fourierprojection(pdegree, theta)
    # k = numbasis(s, tol)
    # U = U[:,1:k]

    return U, s 
end


function cosinprojection(pdegree, theta, theta0)
    
    N = length(theta);
    M = length(theta0);
    nbf = pdegree*M;
    abf = zeros(N, nbf);
    
    for j = 1:M    
        for i = 1:pdegree            
            #abf[:,i+(j-1)*pdegree] = (cos.(theta)*cos(theta0[j]) + sin.(theta)*sin(theta0[j])).^i; 
            abf[:,i+(j-1)*pdegree] = cos.(i*theta .+ theta0[j]); 
        end
    end    

    abfnew, U, s = eigendecomposition(abf, theta)

    return abf, abfnew, U, s
end

function cosinbasisprojection(pdegree, theta0, tol)

    theta = collect(LinRange(0.0, pi, 2000))     
    rbf, rbfnew, U, s = cosinprojection(pdegree, theta, theta0)
    
    return U, s
end


function eigendecomposition(B, theta)

    N = size(B,1)
    A = (1.0/N)*(B'*B);
    s = eigvals(A)
    U = eigvecs(A)
    ind = sortperm(s, rev=true)
    s = s[ind]
    U = U[:,ind]

    Bnew = B*U;
    sij = theta[2:end] - theta[1:end-1];
    for i = 1:size(Bnew,2)
        B2 = Bnew[:,i].*Bnew[:,i];        
        area = 0.5*sum(sij.*(B2[2:end]+B2[1:end-1]));    
        U[:,i] = U[:,i]/sqrt(area);
    end
    Bnew = B*U; 

    return Bnew, U, s 
end

function numbasis(s, tol=1e-4)
    nbf = length(s)
    sumeig = sum(s);
    k = 0
    for i = 1:nbf
        a = sum(s[1:i])/sumeig;
        if a >= (1.0 - tol)
            k = i
            break;
        end
    end
    return k
end

function radialbasis2(rij, rin, rcut, scalefac, pdegree, K)

    epsil = 1e-6;    
    dim, N = size(rij);
    M = length(scalefac);
    nbf = pdegree*M + K
    rbf = zeros(N, nbf);
    drbf = zeros(dim*N, nbf);
    
    rmax = rcut - rin
    for n = 1:N
        dij = sqrt(rij[1,n]*rij[1,n] + rij[2,n]*rij[2,n] + rij[3,n]*rij[3,n])    
        #dr = rij[:,n]/dij;        
        dr1 = rij[1,n]/dij;    
        dr2 = rij[2,n]/dij;    
        dr3 = rij[3,n]/dij;    
        r = dij - rin        
        #fcut,dfcut = cutoffderiv(dij, rin, rcut, epsil) 
        y = r/rmax;    
        fcut = exp(-1.0/sqrt((1.0 - y^3)^2 + epsil))/exp(-1.0);
        dfcut = ((3.0/(rmax*exp(-1.0)))*(y^2)*exp(-1.0/(epsil + (y^3 - 1.0)^2)^(1/2))*(y^3 - 1.0))/(epsil + (y^3 - 1.0)^2)^(3/2);
        #fcut = 0.5*(1.0 .+ cos.(pi*y))
        #dfcut = -(pi*sin.(pi*y))/(2.0*rmax)

        for j = 1:M
            alpha = scalefac[j];    
            if alpha == 0
                alpha = 1e-3;
            end
            x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
            dx = (alpha/rmax)*exp(-(alpha*r/rmax))/(1.0 - exp(-alpha));        

            for i = 1:pdegree    
                a = i*pi;
                rbf[n,i+(j-1)*pdegree] = (sqrt(2.0/(rmax))/i)*fcut*sin(a*x)/r;
                drbfdr = (sqrt(2/rmax)/i)*(dfcut*sin(a*x)/r - fcut*sin(a*x)/(r^2) + a*cos(a*x)*fcut*dx/r);
                #drbf[:,n,i+(j-1)*pdegree] = drbfdr*dr;
                drbf[1+3*(n-1),i+(j-1)*pdegree] = drbfdr*dr1;
                drbf[2+3*(n-1),i+(j-1)*pdegree] = drbfdr*dr2;
                drbf[3+3*(n-1),i+(j-1)*pdegree] = drbfdr*dr3;
            end
        end

        # add radial basis
        for i = 1:K
            p = pdegree*M+i;
            rbf[n,p] = fcut/(dij^i);
            drbfdr = dfcut/(dij^i) - i*fcut/(dij^(i+1));  
            #drbf[:,n,p] = drbfdr*dr;
            drbf[1+3*(n-1),p] = drbfdr*dr1;
            drbf[2+3*(n-1),p] = drbfdr*dr2;
            drbf[3+3*(n-1),p] = drbfdr*dr3;
        end
    end
    
    return rbf, drbf 
end

function onebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)

    dim, N = size(x)
    M = length(pod.species)
    d = ones(M,1)
    for m = 1:M
        ind = findall(t[:] .== m);
        d[m] = 1.0*length(ind)
    end
    dd = zeros(dim,N,M)
    dv = ones(6,M)

    # dim, N = size(x)
    # typei = pod.atomtypes[1];    
    # ilist = findatomtype(Array(1:N), t, typei);      
    # xi = x[:,ilist]
    # ai = ilist
    # ti = t[ilist]    
    # ei = ones(length(ti))
    # fi = zeros(size(xi))        
    # ea, fa = tallysingle(ei, fi, ai);   
    # d = [sum(ea)]         
    # dd = fa
    # va = tallysinglevirial(fi, xi, ai, -1.0)
    # dv = sum(va, dims=2)
    # dd = reshape(dd, (dim,size(fa,2),1))        
    return d, dd, dv 
end

function twobodydescriptors(x, t, a, b, c, pbc, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    S2 = Int32(S*(S+1)/2)
    M = size(pod.Phi,2)
    eatom = zeros(natom,M*S2)
    fatom = zeros(3,natom,M*S2)
    vatom = zeros(6,natom,M*S2)    

    rcut = pod.rcut 
    y, alist, neighlist, pairnumsum, pairnum = fullneighborlist(x, a, b, c, pbc, rcut);

    Nij = length(neighlist)        
    if Nij > 0
        pairlist = Int32.(neighlist[:].-1)
        alist = Int32.(alist.-1)
        atomtype = Int32.(t[:])        

        rij = zeros(3,Nij) 
        ai = zeros(Int32,Nij)    
        aj = zeros(Int32,Nij)  
        ti = zeros(Int32,Nij)      
        tj = zeros(Int32,Nij)        
        NeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
            alist, Int32(natom), Int32(dim))
    
        Ngm = length(pod.gamma0)
        nbf = pod.pdegree[1]*Ngm + pod.pdegree[2]
        eij = zeros(Nij, nbf);
        fij = zeros(dim*Nij, nbf);    
        radialbasis(eij, fij, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, 
            Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
        eij = eij*pod.Phi 
        fij = fij*pod.Phi         

        ind = zeros(Int32,S,S)
        k = 0
        for i = 1:S
            for j = i:S
                k = k + 1
                ind[i,j] = k
                ind[j,i] = k
            end
        end
    
        N,M = size(eij)    
        for n = 1:N
            i1 = ai[n]+1
            j1 = aj[n]+1
            typei = ti[n]
            typej = tj[n]
            mij = (ind[typei,typej]-1)*M
            for m = 1:M  
                ms = m + mij 
                eatom[i1,ms] += eij[n,m]
                nn = 3*(n-1)
                fatom[1,i1,ms] += fij[1+nn,m]
                fatom[2,i1,ms] += fij[2+nn,m]
                fatom[3,i1,ms] += fij[3+nn,m]
                fatom[1,j1,ms] -= fij[1+nn,m]
                fatom[2,j1,ms] -= fij[2+nn,m]
                fatom[3,j1,ms] -= fij[3+nn,m]
    
                v0 = rij[1,n]*fij[1+nn,m]        
                v1 = rij[2,n]*fij[2+nn,m]        
                v2 = rij[3,n]*fij[3+nn,m]        
                v3 = rij[2,n]*fij[3+nn,m]        
                v4 = rij[1,n]*fij[3+nn,m]        
                v5 = rij[1,n]*fij[2+nn,m]      
                vatom[1,i1,ms] += v0; 
                vatom[2,i1,ms] += v1;
                vatom[3,i1,ms] += v2; 
                vatom[4,i1,ms] += v3;
                vatom[5,i1,ms] += v4; 
                vatom[6,i1,ms] += v5;        
                vatom[1,j1,ms] += v0; 
                vatom[2,j1,ms] += v1;
                vatom[3,j1,ms] += v2; 
                vatom[4,j1,ms] += v3;
                vatom[5,j1,ms] += v4; 
                vatom[6,j1,ms] += v5;        
            end
        end    
    end

    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 
    dv = (-1.0/2.0)*reshape(sum(vatom, dims=2), (6,M*S2))

    return d, dd, dv 
end

function twobodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    S2 = Int32(S*(S+1)/2)
    M = size(pod.Phi,2)
    eatom = zeros(natom,M*S2)
    fatom = zeros(3,natom,M*S2)
    vatom = zeros(6,natom,M*S2)    

    ind = zeros(Int32,S,S)
    k = 0
    for i = 1:S
        for j = i:S
            k = k + 1
            ind[i,j] = k
            ind[j,i] = k
        end
    end

    rcut = pod.rcut 
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:natom))

    pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut)
    rij, ai, aj, ti, tj = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);            

    eij, fij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    eij = eij*pod.Phi 
    fij = fij*pod.Phi     
    
    N = size(eij,1)    
    for n = 1:N
        i1 = ai[n]
        j1 = aj[n]
        typei = ti[n]
        typej = tj[n]
        mij = (ind[typei,typej]-1)*M
        for m = 1:M  
            ms = m + mij 
            eatom[i1,ms] += eij[n,m]
            nn = 3*(n-1)
            fatom[1,i1,ms] += fij[1+nn,m]
            fatom[2,i1,ms] += fij[2+nn,m]
            fatom[3,i1,ms] += fij[3+nn,m]
            fatom[1,j1,ms] -= fij[1+nn,m]
            fatom[2,j1,ms] -= fij[2+nn,m]
            fatom[3,j1,ms] -= fij[3+nn,m]

            v0 = rij[1,n]*fij[1+nn,m]        
            v1 = rij[2,n]*fij[2+nn,m]        
            v2 = rij[3,n]*fij[3+nn,m]        
            v3 = rij[2,n]*fij[3+nn,m]        
            v4 = rij[1,n]*fij[3+nn,m]        
            v5 = rij[1,n]*fij[2+nn,m]      
            vatom[1,i1,ms] += v0; 
            vatom[2,i1,ms] += v1;
            vatom[3,i1,ms] += v2; 
            vatom[4,i1,ms] += v3;
            vatom[5,i1,ms] += v4; 
            vatom[6,i1,ms] += v5;        
            vatom[1,j1,ms] += v0; 
            vatom[2,j1,ms] += v1;
            vatom[3,j1,ms] += v2; 
            vatom[4,j1,ms] += v3;
            vatom[5,j1,ms] += v4; 
            vatom[6,j1,ms] += v5;        
        end
    end
    
    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 
    dv = (-1.0/2.0)*reshape(sum(vatom, dims=2), (6,M*S2))

    return d, dd, dv 
end

function radialdescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, N = size(x)
    rcut = pod.rcut 
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);

    typei = pod.atomtypes[1];
    typej = pod.atomtypes[2];            
    ilist = findatomtype(Array(1:N), t, typei);      
    pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, typej, neighlist, neighnum, rcut*rcut);                                                
    # ilist = Int32.(Array(1:N));                        
    # pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);

    rij, ai, aj = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);            
    
    eij, fij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    eij = eij*pod.Phi 
    fij = fij*pod.Phi 
    M = size(pod.Phi,2)
    
    N = size(eij,1)    
    natom = size(x,2)
    eatom = zeros(natom,M)
    fatom = zeros(3,natom,M)
    vatom = zeros(6,natom,M)    
    for n = 1:N
        i1 = ai[n]
        j1 = aj[n]
        for m = 1:M            
            eatom[i1,m] += eij[n,m]
            nn = 3*(n-1)
            fatom[1,i1,m] += fij[1+nn,m]
            fatom[2,i1,m] += fij[2+nn,m]
            fatom[3,i1,m] += fij[3+nn,m]
            fatom[1,j1,m] -= fij[1+nn,m]
            fatom[2,j1,m] -= fij[2+nn,m]
            fatom[3,j1,m] -= fij[3+nn,m]

            v0 = rij[1,n]*fij[1+nn,m]        
            v1 = rij[2,n]*fij[2+nn,m]        
            v2 = rij[3,n]*fij[3+nn,m]        
            v3 = rij[2,n]*fij[3+nn,m]        
            v4 = rij[1,n]*fij[3+nn,m]        
            v5 = rij[1,n]*fij[2+nn,m]      
            vatom[1,i1,m] += v0; 
            vatom[2,i1,m] += v1;
            vatom[3,i1,m] += v2; 
            vatom[4,i1,m] += v3;
            vatom[5,i1,m] += v4; 
            vatom[6,i1,m] += v5;        
            vatom[1,j1,m] += v0; 
            vatom[2,j1,m] += v1;
            vatom[3,j1,m] += v2; 
            vatom[4,j1,m] += v3;
            vatom[5,j1,m] += v4; 
            vatom[6,j1,m] += v5;        
        end
    end
    
    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 
    dv = (-1.0/2.0)*reshape(sum(vatom, dims=2), (6,M))

    # d2 = sum(eatom.^2, dims=1)
    # d2 = d2[:]
    # dd2 = reshape(2*eatom, (1, natom, M)) .* fatom 
    # dv2 = (-1.0/2.0)*reshape(sum(vatom, dims=2), (6,M))
    # d = [d[:]; d2]
    # dd = cat(dd, dd2, dims=3)
    # dv = [dv dv2]

    return d, dd, dv 
end

function cosinbasis!(abf, dabf, xij, xik, nabf, theta0)

    dim, N = size(xij);
    n0 = length(theta0)
    # abf  = zeros(N, n0*nabf)
    # dabf  = zeros(6, N, n0*nabf)        

    for n = 1:N
        xij1 = xij[1,n]
        xij2 = xij[2,n]
        xij3 = xij[3,n]
        xik1 = xik[1,n]
        xik2 = xik[2,n]
        xik3 = xik[3,n]

        xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;
        rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
        riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;    
        rij = sqrt(rijsq); 
        rik = sqrt(riksq); 

        costhe = xdot/(rij*rik);    
        #ind = findall(abs.(costhe) .>= 1.0)
        #costhe[ind] = 1.0*sign.(costhe[ind])
        if abs(costhe) > 1.0
            costhe = 1.0*sign(costhe)
        end 
        xdot = costhe*(rij*rik)

        sinthe = sqrt(1 - costhe*costhe);
        sinthe = max(sinthe, 1e-12)    
        theta = acos(costhe)            
        dtheta = -1.0/sinthe 

        tm1 = (rijsq^(3/2))*(riksq^(1/2));
        tm2 = (rijsq^(1/2))*(riksq^(3/2));
        dct1 = (xik1*rijsq - xij1*xdot)/tm1; 
        dct2 = (xik2*rijsq - xij2*xdot)/tm1;
        dct3 = (xik3*rijsq - xij3*xdot)/tm1;
        dct4 = (xij1*riksq - xik1*xdot)/tm2;
        dct5 = (xij2*riksq - xik2*xdot)/tm2;
        dct6 = (xij3*riksq - xik3*xdot)/tm2;

        for j = 1:n0
            for p = 1:nabf        
                q = p + (j-1)*nabf
                abf[n,q] = cos(p*theta + theta0[j]);                
                tm = -p*sin(p*theta + theta0[j])*dtheta
                dabf[1,n,q] = tm*dct1;
                dabf[2,n,q] = tm*dct2;
                dabf[3,n,q] = tm*dct3;
                dabf[4,n,q] = tm*dct4;
                dabf[5,n,q] = tm*dct5;
                dabf[6,n,q] = tm*dct6;
            end
        end
    end
end

function threebodydescriptors(x, t, a, b, c, pbc, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    nrbf = size(pod.Phi,2)
    nabf = pod.pdegree[3]
    M = nrbf*(nabf+1)  
    S2 = Int32((S+1)*S/2)    

    eatom = zeros(natom,M*S2*S)
    fatom = zeros(3,natom,M*S2*S)

    rcut = pod.rcut 
    y, alist, neighlist, pairnumsum, pairnum = fullneighborlist(x, a, b, c, pbc, rcut);

    Nij = length(neighlist)        
    if Nij > 0
        pairlist = Int32.(neighlist[:].-1)
        alist = Int32.(alist.-1)
        atomtype = Int32.(t[:])        

        ngm = length(pod.gamma0)
        nbf = pod.pdegree[1]*ngm + pod.pdegree[2]
        nelements = length(pod.species)
        Nj = maximum(pairnumsum[2:end]-pairnumsum[1:end-1]);  
        Nijk  = Int32(sum((pairnum.-1).*pairnum/2));

        n = 3*Nij+ (1+dim)*Nij*(nrbf+nbf) + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf + 2*dim*Nijk + 7*Nijk*nabf;
        m = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj + 6*Nijk;
        tmpmem = zeros(Float64,n)
        tmpint = zeros(Int32,m)

        pod3body(eatom, fatom, y, pod.Phi, pod.gamma0, tmpmem, pod.rin, pod.rcut, tmpint, 
            pod.indelem, pairlist, pairnum, pairnumsum, atomtype, alist, Int32.(pod.pdegree), 
            Int32(ngm), Int32(nbf), Int32(nrbf), Int32(nabf), Int32(nelements), Int32(natom), 
            Int32(Nj), Int32(Nij), Int32(Nijk))        
    end            

    dim,N,M = size(fatom)
    d = zeros(M)
    dv = zeros(6,M)
    for m = 1:M
        for n = 1:N
            d[m] += eatom[n,m] 
            vi = (x[:,n])*(fatom[:,n,m]')
            dv[:,m] +=  [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]];
        end
    end

    return d, fatom, dv     
end


function threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc, tm)
        
    dim, natom = size(x)
    S = length(pod.species)
    nrbf = size(pod.Phi,2)
    nabf = pod.pdegree[3]
    M = nrbf*(nabf+1)  
    S2 = Int32((S+1)*S/2)
    t = Int32.(t)

    eatom = zeros(natom,M*S2*S)
    fatom = zeros(3,natom,M*S2*S)

    ind = zeros(Int32,S,S)
    k = 0
    for i = 1:S
        for j = i:S
            k = k + 1
            ind[i,j] = k
            ind[j,i] = k
        end
    end

    #display("start")    
    #@time begin
    rcut = pod.rcut 
    y, alist, neighlist, pairnumsum, pairnum = fullneighborlist(x, a, b, c, pbc, rcut);
    #ilist = Int32.(Array(1:natom))
    #end

    # pairnum = zeros(Int32, natom)
    # pairnumsum = zeros(Int32, natom+1)
    # pairlist = zeros(Int32, length(neighlist))
    # Nij = NeighPairList(pairnum, pairnumsum, pairlist, y, rcut*rcut, Int32.(neighlist.-1), 
    #     Int32.(neighnumsum), Int32(natom), Int32(dim))    

    # display(maximum(abs.(pairnum[:]-neighnum[:]))) 
    # display(maximum(abs.(pairnumsum[:]-neighnumsum[:]))) 
    # display(maximum(abs.(neighlist[:]-pairlist[:].-1)))

    Nij = length(neighlist)        
    if Nij > 0
        pairlist = Int32.(neighlist[:].-1)
        #pairnumsum = Int32.(neighnumsum)
        #pairnum = Int32.(neighnum)
        alist = Int32.(alist.-1)
        atomtype = Int32.(t[:])        
        
        ngm = length(pod.gamma0)
        nbf = pod.pdegree[1]*ngm + pod.pdegree[2]
        nelements = length(pod.species)
        Nj = maximum(pairnumsum[2:end]-pairnumsum[1:end-1]);  
        Nijk  = Int32(sum((pairnum.-1).*pairnum/2));

        n = 3*Nij+ (1+dim)*Nij*(nrbf+nbf) + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf + 2*dim*Nijk + 7*Nijk*nabf;
        m = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj + 6*Nijk;
        tmpmem = zeros(Float64,n)
        tmpint = zeros(Int32,m)

        pod3body(eatom, fatom, y, pod.Phi, pod.gamma0, tmpmem, pod.rin, pod.rcut, tmpint, 
            ind[:], pairlist, pairnum, pairnumsum, atomtype, alist, Int32.(pod.pdegree), 
            Int32(ngm), Int32(nbf), Int32(nrbf), Int32(nabf), Int32(nelements), Int32(natom), 
            Int32(Nj), Int32(Nij), Int32(Nijk))        

        # rij = zeros(3,Nij) 
        # ai = zeros(Int32,Nij)    
        # aj = zeros(Int32,Nij)  
        # ti = zeros(Int32,Nij)      
        # tj = zeros(Int32,Nij)        
        # NeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
        #     alist, Int32(natom), Int32(dim))

        # tripletnum = zeros(Int32, natom)
        # tripletnumsum = zeros(Int32, natom+1)
        # tripletlist = zeros(Int32, Nijk*2)
        # NeighTripletList(tripletlist, tripletnum, tripletnumsum, pairlist, 
        #     pairnum, pairnumsum, Int32(natom))        

        # display("start")    
        # display([a b c])
        # display(rij)
        # display([natom Nj Nij Nijk])
    
        # # #rij0 = tmpmem[1:(3*Nij)]
        # e2ij0 = tmpmem[(3*Nij+1):(3*Nij+Nij*nrbf)];    
        # f2ij0 = tmpmem[(3*Nij+Nij*nrbf+1):(3*Nij+4*Nij*nrbf)];    
        # #e2ij1 = tmpmem[(3*Nij+4*Nij*nrbf+1):(3*Nij+4*Nij*nrbf+Nij*nbf)];    
        # #f2ij1 = tmpmem[(3*Nij+4*Nij*nrbf+Nij*nbf+1):(3*Nij+4*Nij*nrbf+4*Nij*nbf)];  
        # uij0 = tmpmem[(3*Nij+4*Nij*nrbf+1):(3*Nij+4*Nij*nrbf+Nijk*nrbf)];    
        # uik0 = tmpmem[(3*Nij+4*Nij*nrbf+Nijk*nrbf+1):(3*Nij+4*Nij*nrbf+2*Nijk*nrbf)];    
        # wij0 = tmpmem[(3*Nij+4*Nij*nrbf+2*Nijk*nrbf+1):(3*Nij+4*Nij*nrbf+5*Nijk*nrbf)];    
        # wik0 = tmpmem[(3*Nij+4*Nij*nrbf+5*Nijk*nrbf+1):(3*Nij+4*Nij*nrbf+8*Nijk*nrbf)];    
        # xij0 = tmpmem[(3*Nij+4*Nij*nrbf+8*Nijk*nrbf+4*Nj*nrbf+1):(3*Nij+4*Nij*nrbf+8*Nijk*nrbf+3*Nijk+4*Nj*nrbf)];    
        # xik0 = tmpmem[(3*Nij+4*Nij*nrbf+8*Nijk*nrbf+3*Nijk+4*Nj*nrbf+1):(3*Nij+4*Nij*nrbf+8*Nijk*nrbf+6*Nijk+4*Nj*nrbf)];    

        # tripletnum0 = tmpint[(4*Nij+1):(4*Nij+natom)]
        # tripletnumsum0 = tmpint[(4*Nij+natom+1):(4*Nij+2*natom+1)]
        # tripletlist0 = tmpint[(4*Nij+2*natom+2):(4*Nij+2*natom+1+2*Nijk)]
    
        # Ngm = length(pod.gamma0)
        # nbf = pod.pdegree[1]*Ngm + pod.pdegree[2]
        # e2ij = zeros(Nij, nbf);
        # f2ij = zeros(dim*Nij, nbf);    
        # radialbasis(e2ij, f2ij, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, 
        #     Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
        
        # #e2ijt = 1.0*e2ij 
        # #f2ijt = 1.0*f2ij 
        
        # e2ij = e2ij*pod.Phi 
        # f2ij =  f2ij*pod.Phi 

        # M = size(e2ij,2)
        # uij = zeros(Nijk,M);
        # uik = zeros(Nijk,M);      
        # wij = zeros(3,Nijk,M);
        # wik = zeros(3,Nijk,M);              
        # m = Int32((Nj-1)*Nj/2);      
        # indj = zeros(Int32, m)
        # indk = zeros(Int32, m)
        # ei = zeros(Nj, M)
        # f1 = zeros(Nj, M)
        # f2 = zeros(Nj, M)
        # f3 = zeros(Nj, M)
        # makejk(uij, uik, wij, wik, e2ij, f2ij, ei, f1, f2, f3, pairnumsum, tripletnumsum, 
        #     indj, indk, Int32(Nj), Int32(M), Int32(natom), Int32(Nij), Int32(Nijk));

        # xij = zeros(3,Nijk) 
        # xik = zeros(3,Nijk) 
        # ai = zeros(Int32,Nijk)    
        # aj = zeros(Int32,Nijk)  
        # ak = zeros(Int32,Nijk)  
        # ti = zeros(Int32,Nijk)      
        # tj = zeros(Int32,Nijk)   
        # tk = zeros(Int32,Nijk)       
        # NeighTriplets(xij, xik, y, ai, aj, ak, ti, tj, tk, tripletlist, tripletnumsum, 
        #     alist, atomtype,  Int32(natom), Int32(dim))                    

        # display(maximum(abs.(e2ij0[:]-e2ij[:]))) 
        # display(maximum(abs.(f2ij0[:]-f2ij[:])))   
        # #display(maximum(abs.(e2ij1[:]-e2ijt[:]))) 
        # #display(maximum(abs.(f2ij1[:]-f2ijt[:])))                 
        # display(maximum(abs.(uij0[:]-uij[:]))) 
        # display(maximum(abs.(uik0[:]-uik[:]))) 
        # display(maximum(abs.(wij0[:]-wij[:]))) 
        # display(maximum(abs.(wik0[:]-wik[:]))) 
        # display(maximum(abs.(xij0[:]-xij[:]))) 
        # display(maximum(abs.(xik0[:]-xik[:]))) 
        # display(maximum(abs.(tripletnum0[:]-tripletnum[:])))
        # display(maximum(abs.(tripletnumsum0[:]-tripletnumsum[:]))) 
        # display(maximum(abs.(tripletlist0[:]-tripletlist[:])))
        # display("here")
    end    
        

    # Nj = maximum(pairnumsum[2:end]-pairnumsum[1:end-1]);  
    # Nijk  = Int32(sum((pairnum.-1).*pairnum/2));
    # if Nij > 0        
    #     alist = Int32.(alist.-1)
    #     atomtype = Int32.(t[:])        

    #     rij = zeros(3,Nij) 
    #     ai = zeros(Int32,Nij)    
    #     aj = zeros(Int32,Nij)  
    #     ti = zeros(Int32,Nij)      
    #     tj = zeros(Int32,Nij)        
    #     NeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
    #         alist, Int32(natom), Int32(dim))

    #     tripletnum = zeros(Int32, natom)
    #     tripletnumsum = zeros(Int32, natom+1)
    #     tripletlist = zeros(Int32, Nijk*2)
    #     NeighTripletList(tripletlist, tripletnum, tripletnumsum, pairlist, 
    #         pairnum, pairnumsum, Int32(natom))        
    
    #     Ngm = length(pod.gamma0)
    #     nbf = pod.pdegree[1]*Ngm + pod.pdegree[2]
    #     e2ij = zeros(Nij, nbf);
    #     f2ij = zeros(dim*Nij, nbf);    
    #     radialbasis(e2ij, f2ij, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, 
    #         Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
        
    #     e2ij = e2ij*pod.Phi 
    #     f2ij =  f2ij*pod.Phi 
                
    #     M = size(e2ij,2)
    #     uij = zeros(Nijk,M);
    #     uik = zeros(Nijk,M);      
    #     wij = zeros(3,Nijk,M);
    #     wik = zeros(3,Nijk,M);              
    #     m = Int32((Nj-1)*Nj/2);      
    #     indj = zeros(Int32, m)
    #     indk = zeros(Int32, m)
    #     ei = zeros(Nj, M)
    #     f1 = zeros(Nj, M)
    #     f2 = zeros(Nj, M)
    #     f3 = zeros(Nj, M)
    #     makejk(uij, uik, wij, wik, e2ij, f2ij, ei, f1, f2, f3, pairnumsum, tripletnumsum, 
    #         indj, indk, Int32(Nj), Int32(M), Int32(natom), Int32(Nij), Int32(Nijk));
        
    #     xij = zeros(3,Nijk) 
    #     xik = zeros(3,Nijk) 
    #     ai = zeros(Int32,Nijk)    
    #     aj = zeros(Int32,Nijk)  
    #     ak = zeros(Int32,Nijk)  
    #     ti = zeros(Int32,Nijk)      
    #     tj = zeros(Int32,Nijk)   
    #     tk = zeros(Int32,Nijk)       
    #     NeighTriplets(xij, xik, y, ai, aj, ak, ti, tj, tk, tripletlist, tripletnumsum, 
    #         alist, atomtype,  Int32(natom), Int32(dim))                    
    
    #     uijk  = zeros(Nijk, nabf)
    #     wijk  = zeros(6, Nijk, nabf)        
    #     cosinbasis(uijk, wijk, xij, xik, Int32(nabf), Int32(Nijk))
    
    #     podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #         ti, tj, tk, ind[:], Int32(nrbf), Int32(nabf), Int32(natom), Int32(Nijk), Int32(S));        
    # end

    # pairlist, pairnumsum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut)
    # rij,~ = neighpairs(y, [], pairlist, pairnumsum, t, ilist, alist);            

    #pairnum = pairnumsum[2:end] - pairnumsum[1:end-1]; 
    # Nijk  = Int32(sum((pairnum.-1).*pairnum/2));
    # tripletnum = zeros(Int32, natom)
    # tripletnumsum = zeros(Int32, natom+1)
    # tripletlist = zeros(Int32, Nijk*2)
    # if Nijk>0
    #     NeighTripletList(tripletlist, tripletnum, tripletnumsum, Int32.(pairlist), 
    #         pairnum, pairnumsum, Int32(natom))        
    # end

    # xij = zeros(3,Nijk) 
    # xik = zeros(3,Nijk) 
    # ai = zeros(Int32,Nijk)    
    # aj = zeros(Int32,Nijk)  
    # ak = zeros(Int32,Nijk)  
    # ti = zeros(Int32,Nijk)      
    # tj = zeros(Int32,Nijk)   
    # tk = zeros(Int32,Nijk)       
    # if Nijk>0
    #     NeighTriplets(xij, xik, y, ai, aj, ak, ti, tj, tk, tripletlist, tripletnumsum, 
    #         Int32.(alist.-1), Int32.(t[:]),  Int32(natom), Int32(dim))
    # end

    # xij2 = zeros(3,Nijk) 
    # xik2 = zeros(3,Nijk) 
    # ai2 = zeros(Int32,Nijk)    
    # aj2 = zeros(Int32,Nijk)  
    # ak2 = zeros(Int32,Nijk)  
    # ti2 = zeros(Int32,Nijk)      
    # tj2 = zeros(Int32,Nijk)   
    # tk2 = zeros(Int32,Nijk)    
    # NeighTriplets(xij2, xik2, y, ai2, aj2, ak2, ti2, tj2, tk2, tripletlist, tripletnumsum, 
    #     Int32.(alist.-1), Int32.(t[:]),  Int32(natom), Int32(dim))

    # tripletlist = reshape(tripletlist, (Nijk,2))
    # xij, xik, ai, aj, ak, ti, tj, tk = neightriplets(y, [], tripletlist.+1, tripletnumsum, t, ilist, alist);
    # display(xij2)
    # display(xij)
    # display(tripletnumsum)
    # display(tripletlist)
    # display(Nijk)
    # display(pairnum)
    # display(pairlist)
    # display(natom)
    # if maximum(abs.(xij2[:]-xij[:])) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(xij2[:]-xij[:])) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(xik2[:]-xik[:])) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(ai2[:]-ai[:].+1)) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(aj2[:]-aj[:].+1)) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(ak2[:]-ak[:].+1)) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(ti2[:]-ti[:])) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(tj2[:]-tj[:])) > 1e-10
    #     error("here")     
    # end
    # if maximum(abs.(tk2[:]-tk[:])) > 1e-10
    #     error("here")     
    # end

    # tripletlist2, tripletnumsum2 = neightripletlist(pairlist.+1, pairnumsum, ilist);       
    # display(maximum(abs.(tripletnumsum2[:]-tripletnumsum[:])))    
    # display(maximum(abs.(tripletlist2[:]-tripletlist[:])))            
    # xij2, xik2, ai2, aj2, ak2, ti2, tj2, tk2 = neightriplets(y, [], tripletlist2, tripletnumsum, t, ilist, alist);
    # display(maximum(abs.(xij2[:]-xij[:])))      
    # display(maximum(abs.(xik2[:]-xik[:])))               
    # error("here")
    # tripletlist, tripletnumsum = neightripletlist(pairlist.+1, pairnumsum, ilist);
    # xij, xik, ai, aj, ak, ti, tj, tk = neightriplets(y, [], tripletlist, tripletnumsum, t, ilist, alist);
    
    # Nij = size(rij,2)
    # Nijk = size(xij,2)
    # Ngm = length(pod.gamma0)
    # nbf = pod.pdegree[1]*Ngm + pod.pdegree[2]
    # e2ij = zeros(Nij, nbf);
    # f2ij = zeros(dim*Nij, nbf);    
    # radialbasis(e2ij, f2ij, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, 
    #     Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
    
    # e2ij = e2ij*pod.Phi 
    # nm, nrbf = size(pod.Phi)    
    # f2ij =  f2ij*pod.Phi 
    
    # M = size(e2ij,2)
    # uij = zeros(Nijk,M);
    # uik = zeros(Nijk,M);      
    # wij = zeros(3,Nijk,M);
    # wik = zeros(3,Nijk,M);      
    # Nj = maximum(pairnumsum[2:end]-pairnumsum[1:end-1]);  
    # m = Int32((Nj-1)*Nj/2);      
    # indj = zeros(Int32, m)
    # indk = zeros(Int32, m)
    # ei = zeros(Nj, M)
    # f1 = zeros(Nj, M)
    # f2 = zeros(Nj, M)
    # f3 = zeros(Nj, M)
    # makejk(uij, uik, wij, wik, e2ij, f2ij, ei, f1, f2, f3, pairnumsum, tripletnumsum, 
    #     indj, indk, Int32(Nj), Int32(M), Int32(natom), Int32(Nij), Int32(Nijk));

    # uijk  = zeros(Nijk, nabf)
    # wijk  = zeros(6, Nijk, nabf)        
    # cosinbasis(uijk, wijk, xij, xik, Int32(nabf), Int32(Nijk))

    # ai = ai .+ Int32(1)
    # aj = aj .+ Int32(1)
    # ak = ak .+ Int32(1)
    # podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #     ti, tj, tk, ind[:], Int32(nrbf), Int32(nabf), Int32(natom), Int32(Nijk), Int32(S));

    # #@time begin
    # pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut)
    # rij,~ = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);            
    # #end    

    # pairnum2 = zeros(Int32, natom)
    # pairnumsum2 = zeros(Int32, natom+1)
    # pairlist2 = zeros(Int32, natom*100)
    # Nij = NeighPairList(pairnum2, pairnumsum2, pairlist2, y, rcut*rcut, Int32.(neighlist.-1), 
    #     Int32.(neighnum), Int32(natom), Int32(dim))    
    # display(pairnum2)
    # display([pairnum[:] pairnumsum2[:]])
    # n = length(pairlist)
    # display([pairlist[:] pairlist2[1:n]])
    

    # rij2 = zeros(3,Nij) 
    # ai2 = zeros(Int32,Nij)    
    # aj2 = zeros(Int32,Nij)  
    # ti2 = zeros(Int32,Nij)      
    # tj2 = zeros(Int32,Nij)    
    # NeighPairs(rij2, y, ai2, aj2, ti2, tj2, pairlist2, pairnumsum2, Int32.(t), Int32.(alist.-1), Int32(natom), Int32(dim))
    
    # display(maximum(abs.(rij[:]-rij2[:])))    

    # #@time begin
    # tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
    # xij, xik, ai, aj, ak, ti, tj, tk = neightriplets(y, [], tripletlist, tripletnum, t, ilist, alist);
    # #end
    
    # Nijk  = Int32(sum((pairnum2.-1).*pairnum2/2));
    # tripletnum2 = zeros(Int32, natom)
    # tripletnumsum2 = zeros(Int32, natom+1)
    # tripletlist2 = zeros(Int32, 2, Nijk)
    # xij2 = zeros(3,Nijk) 
    # xik2 = zeros(3,Nijk) 
    # ai2 = zeros(Int32,Nijk)    
    # aj2 = zeros(Int32,Nijk)  
    # ak2 = zeros(Int32,Nijk)  
    # ti2 = zeros(Int32,Nijk)      
    # tj2 = zeros(Int32,Nijk)   
    # tk2 = zeros(Int32,Nijk)    
    # NeighTripletList(tripletlist2, tripletnum2, tripletnumsum2, pairlist2, 
    #     pairnum2, pairnumsum2, Int32(natom))        
    # NeighTriplets(xij2, xik2, y, ai2, aj2, ak2, ti2, tj2, tk2, tripletlist2, tripletnumsum2, 
    #     Int32.(alist.-1), Int32.(t),  Int32(natom), Int32(dim))

    # display([tripletnumsum2[:] tripletnum[:]])
    # display([tripletlist2' tripletlist])
    # display(maximum(abs.(xij[:]-xij2[:])))    
    # display(maximum(abs.(xik[:]-xik2[:])))    
    # display([ak[:] ak2[:]])
    # display([tk[:] tk2[:]])
    # error("here")    

    # #e2ij, f2ij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])

    # dim, Nij = size(rij)
    # Ngm = length(pod.gamma0)
    # nbf = pod.pdegree[1]*Ngm + pod.pdegree[2]
    # e2ij = zeros(Nij, nbf);
    # f2ij = zeros(dim*Nij, nbf);    
    # #@time begin
    # radialbasis(e2ij, f2ij, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
    # #end    
    # # Nij = size(rij,2)
    # # Ngm = length(pod.gamma0)
    # # e2ij2 = 0.0*e2ij
    # # f2ij2 = 0.0*f2ij
    # # radialbasis(e2ij2, f2ij2, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
    # # display(maximum(abs.(e2ij[:]-e2ij2[:])))
    # # display(maximum(abs.(f2ij[:]-f2ij2[:])))
    # # display(e2ij)
    # # display(e2ij2)
    # # error("here")

    # e2ij = e2ij*pod.Phi 
    # nm, nrbf = size(pod.Phi)    
    # f2ij =  f2ij*pod.Phi 

    # #@time begin
    # #uij, uik, wij, wik = makejk(e2ij, f2ij, pairnum, tripletnum, ilist)
    # #end

    # inum = length(ilist);
    # Nijk = sum(tripletnum[ilist.+1]-tripletnum[ilist]);        
    # M = size(e2ij,2)
    # uij = zeros(Nijk,M);
    # uik = zeros(Nijk,M);      
    # wij = zeros(3,Nijk,M);
    # wik = zeros(3,Nijk,M);      
    # Nj = maximum(pairnum[ilist.+1]-pairnum[ilist]);  
    # m = Int32((Nj-1)*Nj/2);      
    # indj = zeros(Int32, m)
    # indk = zeros(Int32, m)
    # ei = zeros(Nj, M)
    # f1 = zeros(Nj, M)
    # f2 = zeros(Nj, M)
    # f3 = zeros(Nj, M)
    # #@time begin
    # makejk(uij, uik, wij, wik, e2ij, f2ij, ei, f1, f2, f3, pairnum, tripletnum, 
    #     indj, indk, Int32(Nj), Int32(M), Int32(inum), Int32(Nij), Int32(Nijk));
    # #end
    # # display([uij[1:1000,:] uij2[1:1000,:]])
    # # display([uik[:] uik2[:]])
    # # display([wij[:] wij2[:]])
    # # display([wik[:] wik2[:]])
    # # display(maximum(abs.(uij[:]-uij2[:])))
    # # display(maximum(abs.(uik[:]-uik2[:])))
    # # display(maximum(abs.(wij[:]-wij2[:])))
    # # display(maximum(abs.(wik[:]-wik2[:])))

    # N = size(xij,2)
    # uijk  = zeros(N, nabf)
    # wijk  = zeros(6, N, nabf)        
    # # @time begin    
    # # cosinbasis!(uijk, wijk, xij, xik, nabf, 0.0)    
    # # end
    # #@time begin    
    # cosinbasis(uijk, wijk, xij, xik, Int32(nabf), Int32(N))
    # #end
    # # uijk2  = zeros(N, nabf)
    # # wijk2  = zeros(6, N, nabf)        
    # # cosinbasis(uijk2, wijk2, xij, xik, Int32(nabf), Int32(N))
    # # display(maximum(abs.(uijk[:]-uijk2[:])))
    # # display(maximum(abs.(wijk[:]-wijk2[:])))
    # # display(uijk)
    # # display(uijk2)
    # # error("here")
    
    # # ai = ai .- Int32(1)
    # # aj = aj .- Int32(1)
    # # ak = ak .- Int32(1)
    # #@time begin
    # podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #     ti, tj, tk, ind, Int32(nrbf), Int32(nabf), Int32(natom), Int32(N), Int32(S));
    # #end
    # #display("here")

    # 1  3.80000000000000  0.27728167148849  0.89492661290164  1.17540757248006 
    
    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 

    dim,N,M = size(dd)
    dv = zeros(6,M)
    for m = 1:M
        for n = 1:N
           vi = (x[:,n])*(dd[:,n,m]')
           dv[:,m] +=  [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]];
        end
    end

    return d, dd, dv 
end

function threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    nrbf = size(pod.Phi,2)
    nabf = pod.pdegree[3]
    M = nrbf*(nabf+1)  
    S2 = Int32((S+1)*S/2)
    t = Int32.(t)

    eatom = zeros(natom,M*S2*S)
    fatom = zeros(3,natom,M*S2*S)

    ind = zeros(Int32,S,S)
    k = 0
    for i = 1:S
        for j = i:S
            k = k + 1
            ind[i,j] = k
            ind[j,i] = k
        end
    end

    rcut = pod.rcut 
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:natom))

    pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut)
    rij,~ = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);            

     tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
    xij, xik, ai, aj, ak, ti, tj, tk = neightriplets(y, [], tripletlist, tripletnum, t, ilist, alist);

    e2ij, f2ij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])

    e2ij = e2ij*pod.Phi 
    nm, nrbf = size(pod.Phi)    
    f2ij =  f2ij*pod.Phi 

    uij, uik, wij, wik = makejk(e2ij, f2ij, pairnum, tripletnum, ilist)

    N = size(xij,2)
    uijk  = zeros(N, nabf)
    wijk  = zeros(6, N, nabf)        
    cosinbasis!(uijk, wijk, xij, xik, nabf, 0.0)    
    
    fij = zeros(3)
    fik = zeros(3)    
    for n = 1:N
        K = 0
        typei = ti[n]
        typej = tj[n]
        typek = tk[n]
        mijk = (ind[typej,typek]-1)*M + (typei-1)*M*S2
        for m = 1:(nrbf)
            K = K + 1
            eijk = uij[n,m]*uik[n,m];
            fij[1] = wij[1,n,m]*uik[n,m];
            fij[2] = wij[2,n,m]*uik[n,m];
            fij[3] = wij[3,n,m]*uik[n,m];
            fik[1] = uij[n,m]*wik[1,n,m];
            fik[2] = uij[n,m]*wik[2,n,m];
            fik[3] = uij[n,m]*wik[3,n,m];            

            ms = K + mijk
            i1 = ai[n]
            eatom[i1,ms] += eijk
            fatom[1,i1,ms] += fij[1] + fik[1]
            fatom[2,i1,ms] += fij[2] + fik[2]
            fatom[3,i1,ms] += fij[3] + fik[3]
            
            i1 = aj[n]
            fatom[1,i1,ms] -= fij[1]
            fatom[2,i1,ms] -= fij[2]
            fatom[3,i1,ms] -= fij[3]

            i1 = ak[n]
            fatom[1,i1,ms] -= fik[1]   
            fatom[2,i1,ms] -= fik[2]   
            fatom[3,i1,ms] -= fik[3]   

            for p = 1:(nabf)                
                K = K + 1                
                ujk = uij[n,m]*uik[n,m]    
                uk = uik[n,m]*uijk[n,p]
                uj = uij[n,m]*uijk[n,p]
                eijk = ujk*uijk[n,p];
                fij[1] = wij[1,n,m]*uk + ujk*wijk[1,n,p];
                fij[2] = wij[2,n,m]*uk + ujk*wijk[2,n,p];
                fij[3] = wij[3,n,m]*uk + ujk*wijk[3,n,p];
                fik[1] = wik[1,n,m]*uj + ujk*wijk[4,n,p];
                fik[2] = wik[2,n,m]*uj + ujk*wijk[5,n,p];
                fik[3] = wik[3,n,m]*uj + ujk*wijk[6,n,p];
                            
                ms = K + mijk 
                i1 = ai[n]
                eatom[i1,ms] += eijk
                fatom[1,i1,ms] += fij[1] + fik[1]
                fatom[2,i1,ms] += fij[2] + fik[2]
                fatom[3,i1,ms] += fij[3] + fik[3]
                
                i1 = aj[n]
                fatom[1,i1,ms] -= fij[1]
                fatom[2,i1,ms] -= fij[2]
                fatom[3,i1,ms] -= fij[3]
    
                i1 = ak[n]
                fatom[1,i1,ms] -= fik[1]   
                fatom[2,i1,ms] -= fik[2]   
                fatom[3,i1,ms] -= fik[3]   
            end
        end
    end
    # s = 0
    # for s1 = 1:S
    #     typej = s1
    #     for s2 = s1:S
    #         typek = s2 
    #         s += 1
    #         pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej; typek], neighlist, neighnum, rcut*rcut);                                                
    #         tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
    #         xij, xik, ai, aj, ak, ti = neightriplets(y, [], tripletlist, tripletnum, t, ilist, alist);
        
    #         rij,~ = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);                
    #         e2ij, f2ij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    #         e2ij = e2ij*pod.Phi 
    #         nm, nrbf = size(pod.Phi)    
    #         f2ij =  f2ij*pod.Phi #reshape(f2ij*pod.Phi, (3, nn, nrbf)) 
    #         uij, uik, wij, wik = makejk(e2ij, f2ij, pairnum, tripletnum, ilist)
            
    #         N = size(xij,2)
    #         uijk  = zeros(N, nabf)
    #         wijk  = zeros(6, N, nabf)        
    #         cosinbasis!(uijk, wijk, xij, xik, nabf, 0.0)
            
    #         fij = zeros(3)
    #         fik = zeros(3)    
    #         for n = 1:N
    #             K = 0
    #             typei = ti[n]
    #             for m = 1:(nrbf)
    #                 K = K + 1
    #                 eijk = uij[n,m]*uik[n,m];
    #                 fij[1] = wij[1,n,m]*uik[n,m];
    #                 fij[2] = wij[2,n,m]*uik[n,m];
    #                 fij[3] = wij[3,n,m]*uik[n,m];
    #                 fik[1] = uij[n,m]*wik[1,n,m];
    #                 fik[2] = uij[n,m]*wik[2,n,m];
    #                 fik[3] = uij[n,m]*wik[3,n,m];            

    #                 im = K + s*M + (typei-1)*M*S2
    #                 i1 = ai[n]
    #                 eatom[i1,im] += eijk
    #                 fatom[1,i1,im] += fij[1] + fik[1]
    #                 fatom[2,i1,im] += fij[2] + fik[2]
    #                 fatom[3,i1,im] += fij[3] + fik[3]
                    
    #                 i1 = aj[n]
    #                 fatom[1,i1,im] -= fij[1]
    #                 fatom[2,i1,im] -= fij[2]
    #                 fatom[3,i1,im] -= fij[3]

    #                 i1 = ak[n]
    #                 fatom[1,i1,im] -= fik[1]   
    #                 fatom[2,i1,im] -= fik[2]   
    #                 fatom[3,i1,im] -= fik[3]   

    #                 for p = 1:(nabf)                
    #                     K = K + 1                
    #                     ujk = uij[n,m]*uik[n,m]    
    #                     uk = uik[n,m]*uijk[n,p]
    #                     uj = uij[n,m]*uijk[n,p]
    #                     eijk = ujk*uijk[n,p];
    #                     fij[1] = wij[1,n,m]*uk + ujk*wijk[1,n,p];
    #                     fij[2] = wij[2,n,m]*uk + ujk*wijk[2,n,p];
    #                     fij[3] = wij[3,n,m]*uk + ujk*wijk[3,n,p];
    #                     fik[1] = wik[1,n,m]*uj + ujk*wijk[4,n,p];
    #                     fik[2] = wik[2,n,m]*uj + ujk*wijk[5,n,p];
    #                     fik[3] = wik[3,n,m]*uj + ujk*wijk[6,n,p];
                                    
    #                     im = K + s*M + (typei-1)*M*S2
    #                     i1 = ai[n]
    #                     eatom[i1,im] += eijk
    #                     fatom[1,i1,im] += fij[1] + fik[1]
    #                     fatom[2,i1,im] += fij[2] + fik[2]
    #                     fatom[3,i1,im] += fij[3] + fik[3]
                        
    #                     i1 = aj[n]
    #                     fatom[1,i1,im] -= fij[1]
    #                     fatom[2,i1,im] -= fij[2]
    #                     fatom[3,i1,im] -= fij[3]
            
    #                     i1 = ak[n]
    #                     fatom[1,i1,im] -= fik[1]   
    #                     fatom[2,i1,im] -= fik[2]   
    #                     fatom[3,i1,im] -= fik[3]   
    #                 end
    #             end
    #         end
        
    #     end
    # end
    
    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 

    dim,N,M = size(dd)
    dv = zeros(6,M)
    for m = 1:M
        for n = 1:N
           vi = (x[:,n])*(dd[:,n,m]')
           dv[:,m] +=  [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]];
        end
    end

    return d, dd, dv 
end

function angulardescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, N = size(x)
    rcut = pod.rcut 
    #@time begin
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    #end        

    #@time begin
    typei = pod.atomtypes[1];
    typej = pod.atomtypes[2];            
    typek = pod.atomtypes[3];            
    ilist = findatomtype(Array(1:N), t, typei);      
    pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej; typek], neighlist, neighnum, rcut*rcut);                                                
    #end
    # ilist = Array(1:N);                     
    # pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);         

    #@time begin
    tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
    xij, xik, ai, aj, ak = neightriplets(y, [], tripletlist, tripletnum, t, ilist, alist);
    #end

    #@time begin
    rij,~ = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);                
    #end
    #@time begin
    e2ij, f2ij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    #end
    #@time begin
    e2ij = e2ij*pod.Phi 
    nm, nrbf = size(pod.Phi)    
    #nn = size(e2ij,1)    
    f2ij =  f2ij*pod.Phi #reshape(f2ij*pod.Phi, (3, nn, nrbf)) 
    uij, uik, wij, wik = makejk(e2ij, f2ij, pairnum, tripletnum, ilist)
    #end
    
    #@time begin    
    nabf = pod.pdegree[3]
    n0 = length(pod.theta0)     
    N = size(xij,2)
    uijk  = zeros(N, n0*nabf)
    wijk  = zeros(6, N, n0*nabf)        
    cosinbasis!(uijk, wijk, xij, xik, pod.pdegree[3], pod.theta0)
    #end
    
    nabf = size(uijk,2)
    M = nrbf*(nabf+1)  
    N = size(uij,1)    
    fij = zeros(3)
    fik = zeros(3)    
    natom = size(x,2)
    eatom = zeros(natom,M)
    fatom = zeros(3,natom,M)
    vatom = zeros(6,natom,M)    
    # @time begin
    # for n = 1:N
    #     K = 0
    #     for m = 1:(nrbf)
    #         K = K + 1
    #         eijk = uij[n,m]*uik[n,m];
    #         fij[1] = wij[1,n,m]*uik[n,m];
    #         fij[2] = wij[2,n,m]*uik[n,m];
    #         fij[3] = wij[3,n,m]*uik[n,m];
    #         fik[1] = uij[n,m]*wik[1,n,m];
    #         fik[2] = uij[n,m]*wik[2,n,m];
    #         fik[3] = uij[n,m]*wik[3,n,m];            

    #         v0 = xij[1,n]*fij[1] + xik[1,n]*fik[1];
    #         v1 = xij[2,n]*fij[2] + xik[2,n]*fik[2];        
    #         v2 = xij[3,n]*fij[3] + xik[3,n]*fik[3];
    #         v3 = xij[2,n]*fij[3] + xik[2,n]*fik[3];
    #         v4 = xij[1,n]*fij[3] + xik[1,n]*fik[3];            
    #         v5 = xij[1,n]*fij[2] + xik[1,n]*fik[2];                                

    #         i1 = ai[n]
    #         eatom[i1,K] += eijk
    #         fatom[1,i1,K] += fij[1] + fik[1]
    #         fatom[2,i1,K] += fij[2] + fik[2]
    #         fatom[3,i1,K] += fij[3] + fik[3]
    #         vatom[1,i1,K] += v0; 
    #         vatom[2,i1,K] += v1;
    #         vatom[3,i1,K] += v2; 
    #         vatom[4,i1,K] += v3;
    #         vatom[5,i1,K] += v4; 
    #         vatom[6,i1,K] += v5;        
            
    #         i1 = aj[n]
    #         fatom[1,i1,K] -= fij[1]
    #         fatom[2,i1,K] -= fij[2]
    #         fatom[3,i1,K] -= fij[3]
    #         vatom[1,i1,K] += v0; 
    #         vatom[2,i1,K] += v1;
    #         vatom[3,i1,K] += v2; 
    #         vatom[4,i1,K] += v3;
    #         vatom[5,i1,K] += v4; 
    #         vatom[6,i1,K] += v5;        

    #         i1 = ak[n]
    #         fatom[1,i1,K] -= fik[1]   
    #         fatom[2,i1,K] -= fik[2]   
    #         fatom[3,i1,K] -= fik[3]   
    #         vatom[1,i1,K] += v0; 
    #         vatom[2,i1,K] += v1;
    #         vatom[3,i1,K] += v2; 
    #         vatom[4,i1,K] += v3;
    #         vatom[5,i1,K] += v4; 
    #         vatom[6,i1,K] += v5;        

    #         for p = 1:(nabf)                
    #             K = K + 1                
    #             # eijk = uij[n,m]*uik[n,m]*uijk[n,p];
    #             # fij[1] = wij[1,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[1,n,p];
    #             # fij[2] = wij[2,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[2,n,p];
    #             # fij[3] = wij[3,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[3,n,p];
    #             # fik[1] = uij[n,m]*wik[1,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[4,n,p];
    #             # fik[2] = uij[n,m]*wik[2,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[5,n,p];
    #             # fik[3] = uij[n,m]*wik[3,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[6,n,p];
    #             ujk = uij[n,m]*uik[n,m]    
    #             uk = uik[n,m]*uijk[n,p]
    #             uj = uij[n,m]*uijk[n,p]
    #             eijk = ujk*uijk[n,p];
    #             fij[1] = wij[1,n,m]*uk + ujk*wijk[1,n,p];
    #             fij[2] = wij[2,n,m]*uk + ujk*wijk[2,n,p];
    #             fij[3] = wij[3,n,m]*uk + ujk*wijk[3,n,p];
    #             fik[1] = wik[1,n,m]*uj + ujk*wijk[4,n,p];
    #             fik[2] = wik[2,n,m]*uj + ujk*wijk[5,n,p];
    #             fik[3] = wik[3,n,m]*uj + ujk*wijk[6,n,p];
                
    #             v0 = xij[1,n]*fij[1] + xik[1,n]*fik[1];
    #             v1 = xij[2,n]*fij[2] + xik[2,n]*fik[2];        
    #             v2 = xij[3,n]*fij[3] + xik[3,n]*fik[3];
    #             v3 = xij[2,n]*fij[3] + xik[2,n]*fik[3];
    #             v4 = xij[1,n]*fij[3] + xik[1,n]*fik[3];            
    #             v5 = xij[1,n]*fij[2] + xik[1,n]*fik[2];                                
    
    #             i1 = ai[n]
    #             eatom[i1,K] += eijk
    #             fatom[1,i1,K] += fij[1] + fik[1]
    #             fatom[2,i1,K] += fij[2] + fik[2]
    #             fatom[3,i1,K] += fij[3] + fik[3]
    #             vatom[1,i1,K] += v0; 
    #             vatom[2,i1,K] += v1;
    #             vatom[3,i1,K] += v2; 
    #             vatom[4,i1,K] += v3;
    #             vatom[5,i1,K] += v4; 
    #             vatom[6,i1,K] += v5;        
                
    #             i1 = aj[n]
    #             fatom[1,i1,K] -= fij[1]
    #             fatom[2,i1,K] -= fij[2]
    #             fatom[3,i1,K] -= fij[3]
    #             vatom[1,i1,K] += v0; 
    #             vatom[2,i1,K] += v1;
    #             vatom[3,i1,K] += v2; 
    #             vatom[4,i1,K] += v3;
    #             vatom[5,i1,K] += v4; 
    #             vatom[6,i1,K] += v5;        
    
    #             i1 = ak[n]
    #             fatom[1,i1,K] -= fik[1]   
    #             fatom[2,i1,K] -= fik[2]   
    #             fatom[3,i1,K] -= fik[3]   
    #             vatom[1,i1,K] += v0; 
    #             vatom[2,i1,K] += v1;
    #             vatom[3,i1,K] += v2; 
    #             vatom[4,i1,K] += v3;
    #             vatom[5,i1,K] += v4; 
    #             vatom[6,i1,K] += v5;                            
    #         end
    #     end
    # end
    # end

    ai = Int32.(ai.-1)
    aj = Int32.(aj.-1)
    ak = Int32.(ak.-1)    
    podtally3(eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, 
            ai, aj, ak, Int32(nrbf), Int32(nabf), Int32(natom), Int32(N));
        
    # eatom2 = zeros(natom,M)
    # fatom2 = zeros(3,natom,M)
    # vatom2 = zeros(6,natom,M)    
    # ai = Int32.(ai.-1)
    # aj = Int32.(aj.-1)
    # ak = Int32.(ak.-1)
    # @time begin
    # podtally3(eatom2, fatom2, vatom2, xij, xik, uij, uik, uijk, wij, wik, wijk, 
    #         ai, aj, ak, Int32(nrbf), Int32(nabf), Int32(natom), Int32(N));
    # end
    # display(maximum(abs.(eatom[:]-eatom2[:])))
    # display(maximum(abs.(fatom[:]-fatom2[:])))
    # display(maximum(abs.(vatom[:]-vatom2[:])))
    #podtally3(eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #     nrbf, nabf, natom, N);
    
    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 
    #dv = (-1.0/3.0)*reshape(sum(vatom, dims=2), (6,M))

    dim,N,M = size(dd)
    dv = zeros(6,M)
    for m = 1:M
        for n = 1:N
           vi = (x[:,n])*(dd[:,n,m]')
           dv[:,m] +=  [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]];
        end
    end

    # d2 = sum(eatom.^2, dims=1)
    # d2 = d2[:]
    # dd2 = reshape(2*eatom, (1, natom, M)) .* fatom 
    # dv2 = (-1.0/3.0)*reshape(sum(vatom, dims=2), (6,M))
    # d = [d[:]; d2]
    # dd = cat(dd, dd2, dims=3)
    # dv = [dv dv2]

    return d, dd, dv 
end

function radialbasis(pod::PODdesc, rij::Matrix{Float64})
    eij, fij = radialbasis2(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    eij = eij*pod.Phi 
    fij = fij*pod.Phi 
    return eij, fij 
end

function radialsphericalharmonicbasis(pod::PODdesc, rij::Matrix{Float64})

    g, dg = radialbasis(pod, rij)        
    Ylm, dYlm = cYlm_de!(pod.sh, rij)

    N, K = size(g)
    dg = reshape(dg, (3,N,K));

    M = size(Ylm,2)
    phi = zeros(Complex{Float64},N,M,K);
    dphi = zeros(Complex{Float64},3,N,M,K);
    for k = 1:K
        for m = 1:M    
            for n = 1:N
                phi[n,m,k] = Ylm[n,m]*g[n,k]
                dphi[1,n,m,k] = dYlm[1,n,m]*g[n,k] + Ylm[n,m]*dg[1,n,k]
                dphi[2,n,m,k] = dYlm[2,n,m]*g[n,k] + Ylm[n,m]*dg[2,n,k]
                dphi[3,n,m,k] = dYlm[3,n,m]*g[n,k] + Ylm[n,m]*dg[3,n,k]
            end
        end
    end

    return phi, dphi 
    #return real(phi), imag(phi), real(dphi), imag(dphi)  
end

function radialsphericalharmonicbasis_sum(phi, ai, N)

    Nij, M, K = size(phi)
    A = zeros(Complex{Float64}, N, M, K)
    for k = 1:K     
        for m = 1:M
            for n = 1:Nij 
                ii = ai[n]
                A[ii,m,k] += phi[n,m,k]
            end
        end
    end              

    return A
end

function powerspectrum(Ar, Ai, dAr, dAi, ai, aj, sh::SHBasis)

    N, M, K = size(Ar)
    K2 = Int32(K*(K+1)/2)
    L = Int32(sqrt(M) - 1)
    ps = zeros(N,L+1,K2);
    kk = 0
    for k = 1:K     
        for k1 = k:K           
            kk = kk + 1 
            for n = 1:N
                for l = 0:L
                    for i = 1:(2*l+1)                    
                        ps[n,l+1,kk] += Ar[n,l*l + i,k]*Ar[n,l*l + i,k1] + Ai[n,l*l + i,k]*Ai[n,l*l + i,k1]
                    end
                end
            end
        end
    end

    if (dAr === nothing) | (dAi === nothing) 
        return ps 
    end

    dim, Nij, M, K = size(dAr)
    dps = zeros(dim,N,L+1,K2);        
    kk = 0
    for k = 1:K     
        for k1 = k:K         
            kk = kk + 1    
            for l = 0:L
                for n = 1:Nij
                    ii = ai[n]
                    ij = aj[n]
                    for d = 1:dim
                        dtmp = 0.0
                        for i = 1:(2*l+1)                    
                            dtmp += dAr[d,n,l*l + i,k]*Ar[ii,l*l + i,k1] + 
                                    Ar[ii,l*l + i,k]*dAr[d,n,l*l + i,k1] + 
                                    dAi[d,n,l*l + i,k]*Ai[ii,l*l + i,k1] +
                                    Ai[ii,l*l + i,k]*dAi[d,n,l*l + i,k1] 
                        end
                        dps[d,ii,l+1,kk] += dtmp
                        dps[d,ij,l+1,kk] -= dtmp 
                    end
                end
            end
        end
    end

    return ps, dps   
end

function bispectrum(Ar, Ai, dAr, dAi, ai, aj, pairnum, sh::SHBasis)

    cg = sh.cg
    indl = sh.indl 
    indm = sh.indm
    rowm = sh.rowm 
    
    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K);    
    #@time begin 
    for k = 1:K     
        for j = 1:J
            l2 = indl[j,1];
            l1 = indl[j,2];
            l  = indl[j,3];                                     
            nm = rowm[j+1]-rowm[j];       
            for n = 1:N 
                tmp = 0.0;
                for i = 1:nm
                    q = rowm[j]+i
                    c = cg[q];
                    m2 = indm[q,1];
                    m1 = indm[q,2];
                    m  = indm[q,3];
                    a1 = Ar[n,m + l + (l*l) + 1,k]
                    b1 = Ai[n,m + l + (l*l) + 1,k]
                    a2 = Ar[n,m1 + l1 + (l1*l1) + 1,k]
                    b2 = Ai[n,m1 + l1 + (l1*l1) + 1,k]
                    a3 = Ar[n,m2 + l2 + (l2*l2) + 1,k]
                    b3 = Ai[n,m2 + l2 + (l2*l2) + 1,k]
                    tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                         
                end                                 
                bs[n,j,k] = tmp
            end
        end
    end
    #end

    if (dAr === nothing) | (dAi === nothing) 
        return bs 
    end
    
    dim, Nij, M, K = size(dAr)
    dbsa = zeros(dim,Nij,J,K);        
    Q = size(indm,1);
    bispectrumderiv(dbsa, Ar, Ai, dAr, dAi, cg, pairnum, indl, indm, rowm, 
        Int32(K), Int32(J), Int32(Q), Int32(M), Int32(N), Int32(Nij));

    # @time begin 
    # for k = 1:K     
    #     for j = 1:J
    #         l2 = indl[j,1];
    #         l1 = indl[j,2];
    #         l  = indl[j,3];                                     
    #         nm = rowm[j+1]-rowm[j];       
    #         for n = 1:N 
    #             n0 = pairnum[n]
    #             nj = pairnum[n+1]-n0
    #             for i = 1:nm
    #                 q = rowm[j]+i
    #                 c = cg[q];
    #                 m2 = indm[q,1];
    #                 m1 = indm[q,2];
    #                 m  = indm[q,3];
    #                 ml = m + l + (l*l) + 1          
    #                 ml1 = m1 + l1 + (l1*l1) + 1     
    #                 ml2 = m2 + l2 + (l2*l2) + 1
    
    #                 a1 = Ar[n,ml,k]
    #                 b1 = Ai[n,ml,k]
    #                 a2 = Ar[n,ml1,k]
    #                 b2 = Ai[n,ml1,k]
    #                 a3 = Ar[n,ml2,k]
    #                 b3 = Ai[n,ml2,k]

    #                 a2a3 = a2*a3
    #                 a1a3 = a1*a3
    #                 a1a2 = a1*a2
    #                 b1b2 = b1*b2
    #                 b1b3 = b1*b3
    #                 b2b3 = b2*b3
    #                 a1b2 = a1*b2
    #                 a1b3 = a1*b3
    #                 a2b3 = a2*b3
    #                 a2b1 = a2*b1
    #                 a3b1 = a3*b1
    #                 a3b2 = a3*b2

    #                 c1 = c*(a2a3-b2b3)
    #                 c2 = c*(a1a3+b1b3)
    #                 c3 = c*(a1a2+b1b2)
    #                 c4 = c*(a2b3+a3b2)
    #                 c5 = c*(a3b1-a1b3)
    #                 c6 = c*(a2b1-a1b2)
    #                 for ij = 1:nj
    #                     n1 = n0 + ij                           
    #                     d = 1;      
    #                     dbsa[d,n1,j,k] += dAr[d,n1,ml,k]*c1 + dAr[d,n1,ml1,k]*c2 + dAr[d,n1,ml2,k]*c3 + 
    #                                       dAi[d,n1,ml,k]*c4 + dAi[d,n1,ml1,k]*c5 + dAi[d,n1,ml2,k]*c6;       

    #                     d = 2;                        
    #                     dbsa[d,n1,j,k] += dAr[d,n1,ml,k]*c1 + dAr[d,n1,ml1,k]*c2 + dAr[d,n1,ml2,k]*c3 + 
    #                                       dAi[d,n1,ml,k]*c4 + dAi[d,n1,ml1,k]*c5 + dAi[d,n1,ml2,k]*c6;       

    #                     d = 3;                        
    #                     dbsa[d,n1,j,k] += dAr[d,n1,ml,k]*c1 + dAr[d,n1,ml1,k]*c2 + dAr[d,n1,ml2,k]*c3 + 
    #                                       dAi[d,n1,ml,k]*c4 + dAi[d,n1,ml1,k]*c5 + dAi[d,n1,ml2,k]*c6;       
    #                 end
    #             end                                 
    #         end
    #     end
    # end
    # end

    # Q = size(indm,1);
    # dbsb = zeros(dim,Nij,J,K);     
    # @time begin    
    # bispectrumderiv(dbsb, Ar, Ai, dAr, dAi, cg, Int32.(pairnum), Int32.(indl), Int32.(indm), Int32.(rowm), 
    #         Int32(K), Int32(J), Int32(Q), Int32(M), Int32(N), Int32(Nij));
    # end
    # display(maximum(abs.(dbsb[:]-dbsa[:])))

    #@time begin 
    dbs = zeros(dim,N,J,K);        
    for n = 1:Nij
        ii = ai[n]
        ij = aj[n]
        for k = 1:K     
            for j = 1:J
                dbs[1,ii,j,k] += dbsa[1,n,j,k]
                dbs[1,ij,j,k] -= dbsa[1,n,j,k]
                dbs[2,ii,j,k] += dbsa[2,n,j,k]
                dbs[2,ij,j,k] -= dbsa[2,n,j,k]                     
                dbs[3,ii,j,k] += dbsa[3,n,j,k]
                dbs[3,ij,j,k] -= dbsa[3,n,j,k]                     
            end
        end
    end
    #end

    return bs, dbs   
end

function bispectrumdescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)

    dim, N = size(x)
    rcut = pod.rcut 

    #@time begin 
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);

    typei = pod.atomtypes[1];
    typej = pod.atomtypes[2];            
    ilist = findatomtype(Array(1:N), t, typei);      
    pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, typej, neighlist, neighnum, rcut*rcut);                                                            
    rij, ai, aj = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);            
    #end

    # display(N)
    # display(pairnum)

    #@time begin 
    phi, dphi = radialsphericalharmonicbasis(pod, rij);
    #end
    #@time begin 
    A = radialsphericalharmonicbasis_sum(phi, ai, N);
    #end
    #@time begin 
    #ps, dps = powerspectrum(real(A), imag(A), real(dphi), imag(dphi), ai, aj, pod.sh);
    bs, dbs = bispectrum(real(A), imag(A), real(dphi), imag(dphi), ai, aj, Int32.(pairnum), pod.sh);    
    #end
    #display("End time");

    # dim,N,J2,K2 = size(dps)
    # dim,N,J,K = size(dbs)
    # bs = cat(reshape(ps, (N, J2*K2)), reshape(bs, (N, J*K)), dims=2)
    # dbs = cat(reshape(dps, (dim, N, J2*K2)), reshape(dbs, (dim, N, J*K)), dims=3)

    dim,N,J,K = size(dbs)
    bs = reshape(bs, (N, J*K))
    dbs = reshape(dbs, (dim, N, J*K))
    d = sum(bs, dims=1)

    M = size(dbs,3)
    dv = zeros(6,M)
    for m = 1:M
        for n = 1:N
           vi = (x[:,n])*(dbs[:,n,m]')
           dv[:,m] +=  [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]];
        end
    end

    return d, dbs, dv
end

function PODdescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
    if pod.nbody == 1
        return onebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod)
    elseif pod.nbody == 2
        # @time begin
        # radialdescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        # end
        #@time begin
        #d, dd, dv = twobodydescriptors(x, t, a, b, c, pbc, pod)
        d, dd, dv = twobodydescriptors(x, t, a, b, c, pbc, rcutmax, pod);
        #end
        return d, dd, dv
    elseif pod.nbody == 3
        # @time begin
        # angulardescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        # end
        # display(t)
        # d1, dd1, dv1 = threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod, 0.0)
        # d2, dd2, dv2 = threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod);
        # display(maximum(abs.(d1[:]-d2[:])))    
        # display(maximum(abs.(dd1[:]-dd2[:])))    
        # display(maximum(abs.(dv1[:]-dv2[:])))    
        # error("here")
        #@time begin        
        #d, dd, dv = threebodydescriptors(x, t, a, b, c, pbc, pod)
        d, dd, dv = threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod);
        #end
        return d, dd, dv
    elseif pod.nbody == 4
        # @time begin
        # bispectrumdescriptors(x, t, a, b, c, pbc, rcutmax, pod)        
        # end
        return bispectrumdescriptors(x, t, a, b, c, pbc, rcutmax, pod)        
    else
        error("nbody is invalid")        
    end    
end

