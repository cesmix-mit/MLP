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
    U::Matrix{Float64} = zeros(1,1) # projection matrix 
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
    U::Matrix{Float64} = zeros(1,1), # angular projection matrix 
    )
        
    if length(nbasis) > 1
        sh = SHBasis(nbasis[2]) 
    else
        sh = SHBasis(1) 
    end
    pod = initPOD(PODdesc(name, nbody, species, atomtypes, pdegree, nbasis, rcut, rin, gamma0, theta0, Phi, Psi, Lambda, Upsilon, sh, U))

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

function eigendecomposition(B)

    N = size(B,1)
    A = (1.0/N)*(B'*B);
    s = eigvals(A)
    U = eigvecs(A)
    ind = sortperm(s, rev=true)
    s = s[ind]
    U = U[:,ind]

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

function onebodyefatom(x, t, a, b, c, pbc, pod::PODdesc)

    dim, N = size(x)
    M = length(pod.species)    
    eatom = zeros(N,M)
    for m = 1:M
        ind = findall(t[:] .== m);
        eatom[ind,m] .= 1.0
    end
    fatom = zeros(dim*N,M)

    return eatom, fatom 
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

function twobodyeatom(x, t, a, b, c, pbc, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    S2 = Int32(S*(S+1)/2)
    M = size(pod.Phi,2)
    eatom = zeros(natom,M*S2)

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
        radialbasis0(eij, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, 
            Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
        eij = eij*pod.Phi 

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
            typei = ti[n]
            typej = tj[n]
            mij = (ind[typei,typej]-1)*M
            for m = 1:M  
                ms = m + mij 
                eatom[i1,ms] += eij[n,m]
            end
        end    
    end

    return eatom  
end

function twobodyefatom(x, t, a, b, c, pbc, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    S2 = Int32(S*(S+1)/2)
    M = size(pod.Phi,2)
    eatom = zeros(natom,M*S2)
    fatom = zeros(3*natom,M*S2)

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
                fatom[1+3*(i1-1),ms] += fij[1+nn,m]
                fatom[2+3*(i1-1),ms] += fij[2+nn,m]
                fatom[3+3*(i1-1),ms] += fij[3+nn,m]
                fatom[1+3*(j1-1),ms] -= fij[1+nn,m]
                fatom[2+3*(j1-1),ms] -= fij[2+nn,m]
                fatom[3+3*(j1-1),ms] -= fij[3+nn,m]
            end
        end    
    end

    return eatom, fatom  
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

        # eij0 = zeros(Nij, nbf);
        # radialbasis0(eij0, rij, pod.gamma0, pod.rin, pod.rcut-pod.rin, 
        #     Int32(pod.pdegree[1]), Int32(pod.pdegree[2]), Int32(Ngm), Int32(Nij))
        # eij0 = eij0*pod.Phi 
        # display(maximum(abs.(eij[:]-eij0[:])))            

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

    # dim, natom, M = size(fatom)
    # d2 = sum(eatom.^2, dims=1)
    # #d2 = sum( 2.0*((1.0)./(1.0 .+ exp.(-eatom)) .- 0.5), dims=1)
    # d2 = d2[:]
    # dd2 = reshape(2*eatom, (1, natom, M)) .* fatom 
    # #dd2 = reshape(2*exp.(-eatom)./((exp.(-eatom) .+ 1).^2), (1, natom, M)) .* fatom 
    # d = [d[:]; d2]
    # dd = cat(fatom, dd2, dims=3)
    # dv = [dv dv]

    # display(eatom)
    # display(fatom[1,:,:])
    # error("here")

    # display(eatom)
    # error("here")

# # activation function 
# f(x) = (1.0)./(1.0 .+ exp(-x)) - 0.5

# # Derivative of activation function 
# df(x) = exp(-x)./((exp(-x) + 1).^2) #f(x).*(1-f(x)) 

    return d, dd, dv, eatom  
end

function twobodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    S2 = Int32(S*(S+1)/2)
    M = size(pod.Phi,2)

    #S2 = S
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
        #mij = (typei-1)*M
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
    # ilist = Int64.(Array(1:N));                        
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

function threebodyeatom(x, t, a, b, c, pbc, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    nrbf = size(pod.Phi,2)
    nabf = pod.pdegree[3]
    M = nrbf*(nabf+1)  
    S2 = Int32((S+1)*S/2)    

    eatom = zeros(natom,M*S2*S)

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

        ind = zeros(Int32,S,S)
        k = 0
        for i = 1:S
            for j = i:S
                k = k + 1
                ind[i,j] = k
                ind[j,i] = k
            end
        end

        pod3body0(eatom, y, pod.Phi, pod.gamma0, tmpmem, pod.rin, pod.rcut, tmpint, 
            ind[:], pairlist, pairnum, pairnumsum, atomtype, alist, Int32.(pod.pdegree), 
            Int32(ngm), Int32(nbf), Int32(nrbf), Int32(nabf), Int32(nelements), Int32(natom), 
            Int32(Nj), Int32(Nij), Int32(Nijk)) 
    end            

    return eatom      
end

function threebodyefatom(x, t, a, b, c, pbc, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    nrbf = size(pod.Phi,2)
    nabf = pod.pdegree[3]
    M = nrbf*(nabf+1)  
    S2 = Int32((S+1)*S/2)    

    eatom = zeros(natom,M*S2*S)
    fatom = zeros(3*natom,M*S2*S)   

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

        ind = zeros(Int32,S,S)
        k = 0
        for i = 1:S
            for j = i:S
                k = k + 1
                ind[i,j] = k
                ind[j,i] = k
            end
        end
         
        pod3body(eatom, fatom, y, pod.Phi, pod.gamma0, tmpmem, pod.rin, pod.rcut, tmpint, 
            ind[:], pairlist, pairnum, pairnumsum, atomtype, alist, Int32.(pod.pdegree), 
            Int32(ngm), Int32(nbf), Int32(nrbf), Int32(nabf), Int32(nelements), Int32(natom), 
            Int32(Nj), Int32(Nij), Int32(Nijk)) 
    end            

    return eatom, fatom       
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

        ind = zeros(Int32,S,S)
        k = 0
        for i = 1:S
            for j = i:S
                k = k + 1
                ind[i,j] = k
                ind[j,i] = k
            end
        end

        pod3body(eatom, fatom, y, pod.Phi, pod.gamma0, tmpmem, pod.rin, pod.rcut, tmpint, 
            ind[:], pairlist, pairnum, pairnumsum, atomtype, alist, Int32.(pod.pdegree), 
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

    return d, fatom, dv, eatom      
end

function threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, natom = size(x)
    S = length(pod.species)
    nrbf = size(pod.Phi,2)
    nabf = pod.pdegree[3]
    M = nrbf*(nabf+1)  
    S2 = Int32((S+1)*S/2)
    t = Int32.(t)

    #S2=1
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
    #cosinbasis(uijk, wijk, xij, xik, Int32(nabf), Int32(N))

    # ai = ai .- Int32(1)
    # aj = aj .- Int32(1)
    # ak = ak .- Int32(1)
    # podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #     ti, tj, tk, ind, Int32(nrbf), Int32(nabf), Int32(natom), Int32(N), Int32(S));

    fij = zeros(3)
    fik = zeros(3)    
    for n = 1:N
        K = 0
        typei = ti[n]
        typej = tj[n]
        typek = tk[n]
        mijk = (ind[typej,typek]-1)*M + (typei-1)*M*S2
        #mijk = (typei-1)*M
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

    # display(eatom)
    # show(stdout, "text/plain", eatom)
    # display(fatom[1,:,:])
    # error("here")    
    
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

    # 1  3.80000000000000  0.27728167148849  0.89492661290164  1.17540757248006 
    # function podtally3b(eatom, fatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #     ti, tj, tk, ind, nrbf::Int32, nabf::Int32, natom::Int32, N::Int32, S::Int32)
    
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
        #return twobodydescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        return twobodydescriptors(x, t, a, b, c, pbc, pod)
    elseif pod.nbody == 3
        # @time begin
        # angulardescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        # end
        #return threebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        return threebodydescriptors(x, t, a, b, c, pbc, pod)
    elseif pod.nbody == 4
        # @time begin
        # bispectrumdescriptors(x, t, a, b, c, pbc, rcutmax, pod)        
        # end
        return bispectrumdescriptors(x, t, a, b, c, pbc, rcutmax, pod)        
    else
        error("nbody is invalid")        
    end    
end

function PODdescriptors(x, t, a, b, c, pbc, pod)

    npod = length(pod)
    nbd1 = 0
    nbd2 = 0
    nbd3 = 0
    d1 = 0.0; d2 = 0.0; d3 = 0.0;
    dd1 = 0.0; dd2 = 0.0; dd3 = 0.0;
    for n = 1:npod
        if pod[n].nbody==1
            d1, dd1 = onebodydescriptors(x, t, a, b, c, pbc, pod[n].rcut, pod[n]);
            nbd1 = 1
        elseif pod[n].nbody==2
            d2, dd2 = twobodydescriptors(x, t, a, b, c, pbc, pod[n])
            # display(pod[n].U)
            # display(size(dd2))
            dim,natom,M2 = size(dd2)
            K2 = size(pod[n].U,2)                        
            d2 = reshape(d2,(1,M2))*pod[n].U 
            d2 = d2[:]
            dd2 = reshape(dd2,(dim*natom,M2))*pod[n].U             
            dd2 = reshape(dd2, (dim,natom,K2))
            nbd2 = 1
        elseif pod[n].nbody==3
            d3, dd3 = threebodydescriptors(x, t, a, b, c, pbc, pod[n])
            # display(pod[n].U)
            # display(size(dd3))
            dim,natom,M3 = size(dd3)
            K3 = size(pod[n].U,2)                        
            d3 = reshape(d3,(1,M3))*pod[n].U 
            d3 = d3[:]
            dd3 = reshape(dd3,(dim*natom,M3))*pod[n].U             
            dd3 = reshape(dd3, (dim,natom,K3))
            nbd3 = 1
        end
    end

    if (nbd1==1) & (nbd2==1) & (nbd3==1)
        d = [d1; d2; d3]
        dd = cat(dd1, dd2, dims=3)
        dd = cat(dd, dd3, dims=3)
    elseif (nbd2==1) & (nbd3==1)
        d = [d2; d3]
        dd = cat(dd2, dd3, dims=3)
    else
        error("POD descriptors are not supported");     
    end

    if (nbd2==1) & (nbd3==1)
        dim = size(dd2,1)
        natom = size(dd2,2)
        M2 = length(d2)
        M3 = length(d3)

        d4 = zeros(M2*M3)
        dd4 = zeros(dim, natom, M2*M3)
        podhybrid23(d4, dd4, d2, d3, dd2, dd3, Int32(M2), Int32(M3), Int32(dim*natom))
        # for m3 = 1:M3
        #     for m2 = 1:M2
        #         d4[m2+M2*(m3-1)] = d2[m2]*d3[m3]                
        #         dd4[:,:,m2+M2*(m3-1)] = d2[m2]*dd3[:,:,m3] + dd2[:,:,m2]*d3[m3]
        #     end
        # end
        d4 = d4[:]        
        d = [d; d4]
        dd = cat(dd, dd4, dims=3)
    end

    dim, N, M = size(dd)
    dv = zeros(6,M)

    # d2, dd2 = podefatom(x, t, a, b, c, pbc, pod)
    # display(maximum(abs.(d[:]-d2[:])))            
    # display(maximum(abs.(dd[:]-dd2[:])))            
    return d, dd, dv
end


function podprojection(data, descriptors, pbc=nothing, a=nothing, b=nothing, c=nothing)

    display("Compute POD projection ...")
    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   
    for i = length(data):-1:1    
        if data[i] === nothing
            deleteat!(data, i)            
        end
    end   

    n = length(data)
    B2 = 0.0
    B3 = 0.0
    K2 = 0
    K3 = 0
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        A2, A3, N2, N3 = podeatom(config, descriptors)        
        B2 = B2 .+ A2
        B3 = B3 .+ A3 
        K2 = K2 + N2
        K3 = K3 + N3
    end    

    if K2 > 0
        B2 = (1.0/K2)*B2    
        U2, s2 = eigenvalues(B2)
        n2 = numbasis(s2, 1e-11);
        U2 = U2[:,1:n2];
        #display(size(U2))
    end

    if K3 > 0
        B3 = (1.0/K3)*B3
        U3, s3 = eigenvalues(B3)    
        n3 = numbasis(s3, 1e-11);    
        U3 = U3[:,1:n3];    
        #display(size(U3))
    end

    for n = 1:length(descriptors)
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
            descriptors[n].U = U2 
        end
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
            descriptors[n].U = U3
        end
    end    
    return descriptors
end

function getconfig(config, natom, ci)

    x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
    t = config.t[:,(natom[ci]+1):natom[ci+1]] # atom types of configuration ci

    if length(config.e) > 0
        e = config.e[ci] 
    else
        e = 0.0;
    end
    if size(config.f,2) > 1
        f = config.f[:,(natom[ci]+1):natom[ci+1]]
    else
        f = 0.0;
    end

    if size(config.a,2) > 1
        a = config.a[:,ci]                        
        b = config.b[:,ci]                        
        c = config.c[:,ci]             
    else
        a = config.a[:,1]                        
        b = config.b[:,1]                        
        c = config.c[:,1]                       
    end
    if size(config.pbc,2) > 1
        pbc = config.pbc[:,ci]           
    else
        pbc = config.pbc[:,1]         
    end

    return x, t, a, b, c, pbc, e, f
end

function podeatom(config, descriptors)

# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

A2 = 0.0
A3 = 0.0
N2 = 0
N3 = 0 
for i = 1:nconfigs
    ci = i
    x, t, a, b, c, pbc = getconfig(config, natom, ci)
    # x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
    # t = config.t[:,(natom[ci]+1):natom[ci+1]] # atom types of configuration ci
    # if length(config.q)>0
    #     q = config.q[:,(natom[ci]+1):natom[ci+1]]
    # else
    #     q = []
    # end
    # if size(config.a,2) > 1
    #     a = config.a[:,ci]                        
    #     b = config.b[:,ci]                        
    #     c = config.c[:,ci]             
    # else
    #     a = config.a[:,1]                        
    #     b = config.b[:,1]                        
    #     c = config.c[:,1]                       
    # end
    # if size(config.pbc,2) > 1
    #     pbc = config.pbc[:,ci]           
    # else
    #     pbc = config.pbc[:,1]         
    # end

    for n = 1:length(descriptors)
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
            eatom = twobodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                            
            A2 = A2 .+ eatom'*eatom 
            N2 = N2 .+ size(eatom,1)
        end
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
            eatom = threebodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                
            A3 = A3 .+ eatom'*eatom 
            N3 = N3 .+ size(eatom,1)
        end
    end
end

return A2, A3, N2, N3

end

function podefatom(x, t, a, b, c, pbc, descriptors)

nbd1 = 0; nbd2 = 0; nbd3 = 0    
eatom1 = 0.0; eatom2 = 0.0; eatom3 = 0.0
fatom1 = 0.0; fatom2 = 0.0; fatom3 = 0.0; fatom23 = 0.0
globd1 = 0.0; globd2 = 0.0; globd3 = 0.0; globd23 = 0.0
for n = 1:length(descriptors)
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==1) 
        eatom1, fatom1 = onebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                                        
        globd1 = sum(eatom1, dims=1)            
        nbd1 = 1
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
        eatom2, fatom2 = twobodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                            
        globd2 = sum(eatom2, dims=1)
        globd2 = globd2*descriptors[n].U 
        fatom2 = fatom2*descriptors[n].U             
        nbd2 = 1
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
        eatom3, fatom3 = threebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                
        globd3 = sum(eatom3, dims=1)
        globd3 = globd3*descriptors[n].U 
        fatom3 = fatom3*descriptors[n].U             
        nbd3 = 1
    end
end    

if (nbd2==1) & (nbd3==1)
    dim, N = size(x)
    M2 = length(globd2)
    M3 = length(globd3)

    globd23 = zeros(1,M2*M3)
    fatom23 = zeros(dim*N, M2*M3)    
    podhybrid23(globd23, fatom23, globd2, globd3, fatom2, fatom3, Int32(M2), Int32(M3), Int32(dim*N))
    # for m3 = 1:M3
    #     for m2 = 1:M2
    #         globd23[m2+M2*(m3-1)] = globd2[m2]*globd3[m3]    
    #         for n = 1:(dim*N)            
    #             fatom23[n,m2+M2*(m3-1)] = globd2[m2]*fatom3[n,m3] + fatom2[n,m2]*globd3[m3]
    #         end
    #     end
    # end
    # d23 = zeros(1,M2*M3)
    # dd23 = zeros(dim*N, M2*M3)    
    # podhybrid23(d23, dd23, globd2, globd3, fatom2, fatom3, Int32(M2), Int32(M3), Int32(dim*N))
    # display(maximum(abs.(d23[:]-globd23[:])))
    # display(maximum(abs.(dd23[:]-fatom23[:])))
end           

if (nbd1 == 1) & (nbd2==1) & (nbd3==1)
    globd = cat(globd1, globd2, dims=2)
    globd = cat(globd, globd3, dims=2)
    fatom = cat(fatom1, fatom2, dims=2)
    fatom = cat(fatom, fatom3, dims=2)
    eatom = cat(eatom1, eatom2, dims=2)
    eatom = cat(eatom, eatom3, dims=2)
elseif (nbd2==1) & (nbd3==1)
    globd = cat(globd2, globd3, dims=2)
    fatom = cat(fatom2, fatom3, dims=2)
    eatom = cat(eatom2, eatom3, dims=2)
elseif (nbd2==1) & (nbd1==1)
    globd = cat(globd1, globd2, dims=2)
    fatom = cat(fatom1, fatom2, dims=2)
    eatom = cat(eatom1, eatom2, dims=2)
elseif (nbd2==1) & (nbd3==1)
    globd = cat(globd2, globd3, dims=2)
    fatom = cat(fatom2, fatom3, dims=2)
    eatom = cat(eatom2, eatom3, dims=2)
elseif (nbd2==1) 
    globd = 1.0*globd2         
    fatom = 1.0*fatom2         
    eatom = 1.0*eatom2         
elseif (nbd3==1) 
    globd = 1.0*globd3         
    fatom = 1.0*fatom3         
    eatom = 1.0*eatom3             
end

if (nbd2==1) & (nbd3==1)
    globd = cat(globd, globd23, dims=2)    
    fatom = cat(fatom, fatom23, dims=2)
end

return globd, fatom, eatom 

end

function podlinearsystem(config, descriptors, normalizeenergy)

# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

matA = 0.0
vecb = 0.0
for i = 1:nconfigs
    ci = i
    normconst = 1.0
    if normalizeenergy==1
        normconst = 1.0/config.natom[ci]
    end    
    x, t, a, b, c, pbc, e, f = getconfig(config, natom, ci)
    globd, fatom = podefatom(x, t, a, b, c, pbc, descriptors)
    
    we2 = (config.we[1]*config.we[1])*(normconst*normconst)
    wf2 = (config.wf[1]*config.wf[1])
    matA = matA .+ (we2*(globd'*globd) + wf2*(fatom'*fatom))    
    vecb = vecb .+ ((we2*e)*(globd') + wf2*(fatom'*f[:]))    
end

return matA, vecb

end

function podleastsq(data, descriptors, Doptions)

    display("Perform POD's linear regression ...")

    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   
    for i = length(data):-1:1    
        if data[i] === nothing
            deleteat!(data, i)            
        end
    end   

    pbc = Doptions.pbc
    normalizeenergy = Doptions.normalizeenergy
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c

    n = length(data)
    matA = 0.0
    vecb = 0.0
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        A1, b1 = podlinearsystem(config, descriptors, normalizeenergy)        
        matA = matA .+ A1
        vecb = vecb .+ b1         
    end    
    
    coeff = matA\vecb 
    return coeff, matA, vecb 
end

function poderrors(config, descriptors, coeff, normalizeenergy)

    # cumalative sum of numbers of atoms 
    nconfigs = length(config.natom)
    natom = [0; cumsum(config.natom[:])];
    
    eerr = zeros(nconfigs)
    fmae = zeros(nconfigs)
    frmse = zeros(nconfigs)
    szbf = zeros(nconfigs)
    for i = 1:nconfigs
        ci = i
        normconst = 1.0
        if normalizeenergy==1
            normconst = 1.0/config.natom[ci]
        end    
        x, t, a, b, c, pbc, e, f = getconfig(config, natom, ci)
        globd, fatom = podefatom(x, t, a, b, c, pbc, descriptors)
        e1 = globd*coeff;
        f1 = fatom*coeff;

        eerr[i] = abs(e1[1] - e)*normconst         

        N = length(f1)
        szbf[i] = N 
        res = (f1 - f[:])    
        fmae[i] = sum(abs.(res))/N     
        ssr = sum(res.^2)
        mse = ssr /N;
        frmse[i] = sqrt(mse) 
    end
    
    emae = sum(abs.(eerr))/nconfigs
    ssr = sum(eerr.^2)
    mse = ssr/nconfigs
    ermse = sqrt(mse) 

    nforces = sum(szbf)
    fmaeave = sum(fmae.*szbf)/nforces
    frmseave =  sqrt(sum((frmse.^2).*szbf)/nforces)    

    return eerr, emae, ermse, fmae, frmse, fmaeave, frmseave, nconfigs, nforces
end
    
function poderroranalysis(data, descriptors, Doptions, coeff)

    display("Calculate POD's errors ...")

    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   
    for i = length(data):-1:1    
        if data[i] === nothing
            deleteat!(data, i)            
        end
    end   

    pbc = Doptions.pbc
    normalizeenergy = Doptions.normalizeenergy
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c

    n = length(data)
    emae = -ones(Float64,n)
    fmae = -ones(Float64,n)
    ermse = -ones(Float64,n)
    frmse = -ones(Float64,n)
    szbe = zeros(Int64, n)
    szbf = zeros(Int64, n)
    eerr = Array{Any}(nothing, n)
    ferr = Array{Any}(nothing, n)
    ferm = Array{Any}(nothing, n)
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end

        e1, e2, e3, e4, e5, e6, e7, nconfigs, nforces = 
                poderrors(config, descriptors, coeff, normalizeenergy)        

        eerr[i] = e1
        ferr[i] = e4 
        ferm[i] = e5 

        szbe[i] = nconfigs
        szbf[i] = nforces
        emae[i] = e2
        ermse[i] = e3 
        fmae[i] = e6
        frmse[i] = e7 
    end    
    
    energyerrors = zeros((n+1),2)    
    forceerrors = zeros((n+1),2)    
    stresserrors = zeros((n+1),2)    

    energyerrors[1,1] = sum(emae.*szbe)/sum(szbe)  # MAE
    energyerrors[1,2] = sqrt(sum((ermse.^2).*szbe)/sum(szbe)) # RMSE 
    energyerrors[2:end,:] = [emae ermse]        

    forceerrors[1,1] = sum(fmae.*szbf)/sum(szbf)
    forceerrors[1,2] =  sqrt(sum((frmse.^2).*szbf)/sum(szbf))    
    forceerrors[2:end,:] = [fmae frmse]   

    return energyerrors, forceerrors, stresserrors, eerr, ferr, ferm  
end

