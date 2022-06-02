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
    hybrid23::Bool = false
    projectiontol::Float64 = 1e-11  
    b2b2::Int64 = 0
    b2b3::Int64 = 0
    b2b4::Int64 = 0
    b3b3::Int64 = 0
    b3b4::Int64 = 0
    b4b4::Int64 = 0
    b2b3b4::Int64 = 0
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
    U::Matrix{Float64} = zeros(1,1), #  projection matrix 
    hybrid23::Bool = false,
    projectiontol::Float64 = 1e-11,
    b2b2::Int64 = 0, 
    b2b3::Int64 = 0, 
    b2b4::Int64 = 0, 
    b3b3::Int64 = 0, 
    b3b4::Int64 = 0, 
    b4b4::Int64 = 0, 
    b2b3b4::Int64 = 0   
    )
        
    if length(nbasis) > 1
        sh = SHBasis(nbasis[2]) 
    else
        sh = SHBasis(1) 
    end
    pod = initPOD(PODdesc(name, nbody, species, atomtypes, pdegree, nbasis, rcut, rin, gamma0, 
                theta0, Phi, Psi, Lambda, Upsilon, sh, U, hybrid23, projectiontol, b2b2, b2b3,
                b2b4, b3b3, b3b4, b4b4, b2b3b4))

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

function cosinbasis!(abf, dabf, xij, xik, nabf, theta0)

    dim, N = size(xij);
    n0 = length(theta0)

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

function radialsphericalharmonicbasis_sum(phi, ai, tj, S, N)

    Nij, M, K = size(phi)
    A = zeros(Complex{Float64}, N, M, K, S)
    for k = 1:K     
        for m = 1:M
            for n = 1:Nij 
                ii = ai[n]
                sj = tj[n]
                A[ii,m,k,sj] += phi[n,m,k]
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

function bispectrum(Ar, Ai, dAr, dAi, ai, aj, pairnumsum, sh::SHBasis)

    cg = sh.cg
    indl = sh.indl 
    indm = sh.indm
    rowm = sh.rowm 
    
    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K);    
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

    if (dAr === nothing) | (dAi === nothing) 
        return bs 
    end
    
    dim, Nij, M, K = size(dAr)
    dbsa = zeros(dim,Nij,J,K);        
    Q = size(indm,1);
    bispectrumderiv(dbsa, Ar, Ai, dAr, dAi, cg, pairnumsum, indl[:], indm[:], rowm, 
        Int32(K), Int32(J), Int32(Q), Int32(M), Int32(N), Int32(Nij));

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

    return bs, dbs   
end

function bispectrumdescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)

    dim, N = size(x)
    rcut = pod.rcut 
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);

    typei = pod.atomtypes[1];
    typej = pod.atomtypes[2];            
    ilist = findatomtype(Array(1:N), t, typei);      
    pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, typej, neighlist, neighnum, rcut*rcut);                                                            
    rij, ai, aj = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);            
  
    phi, dphi = radialsphericalharmonicbasis(pod, rij);
    A = radialsphericalharmonicbasis_sum(phi, ai, N);
    #ps, dps = powerspectrum(real(A), imag(A), real(dphi), imag(dphi), ai, aj, pod.sh);
    bs, dbs = bispectrum(real(A), imag(A), real(dphi), imag(dphi), ai, aj, Int32.(pairnum), pod.sh);    

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

function bispectrum(Ar, Ai, dAr, dAi, t, ai, aj, ti, tj, pairnumsum, pod::PODdesc)

    sh = pod.sh
    S = length(pod.species)

    cg = sh.cg
    indl = sh.indl 
    indm = sh.indm
    rowm = sh.rowm 
    
    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K,S);    
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
                bs[n,j,k,t[n]] = tmp
            end
        end
    end

    bs = reshape(bs, (N, J*K*S)) 
    if (dAr === nothing) | (dAi === nothing) 
        return bs 
    end
    
    dim, Nij, M, K = size(dAr)
    dbsa = zeros(dim,Nij,J,K);        
    Q = size(indm,1);
    bispectrumderiv(dbsa, Ar, Ai, dAr, dAi, cg, pairnumsum, indl[:], indm[:], rowm, 
        Int32(K), Int32(J), Int32(Q), Int32(M), Int32(N), Int32(Nij));

    dbs = zeros(dim,N,J,K,S);        
    for n = 1:Nij
        ii = ai[n]
        ij = aj[n]
        t1 = ti[n]
        for k = 1:K     
            for j = 1:J
                dbs[1,ii,j,k,t1] += dbsa[1,n,j,k]
                dbs[1,ij,j,k,t1] -= dbsa[1,n,j,k]
                dbs[2,ii,j,k,t1] += dbsa[2,n,j,k]
                dbs[2,ij,j,k,t1] -= dbsa[2,n,j,k]                     
                dbs[3,ii,j,k,t1] += dbsa[3,n,j,k]
                dbs[3,ij,j,k,t1] -= dbsa[3,n,j,k]                     
            end
        end
    end
    dbs = reshape(dbs, (dim*N, J*K*S)) 

    return bs, dbs   
end

function bispectrum2(Ar, Ai, dAr, dAi, t, ai, aj, ti, tj, pairnumsum, pod::PODdesc)

    sh = pod.sh
    S = length(pod.species)

    cg = sh.cg
    indl = sh.indl 
    indm = sh.indm
    rowm = sh.rowm 
    
    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K,S,S);    
    for k = 1:K     
        for j = 1:J
            l2 = indl[j,1];
            l1 = indl[j,2];
            l  = indl[j,3];                                     
            nm = rowm[j+1]-rowm[j];       
            for n = 1:N 
                for sj = 1:S
                    tmp = 0.0;                
                    for i = 1:nm
                        q = rowm[j]+i
                        c = cg[q];
                        m2 = indm[q,1];
                        m1 = indm[q,2];
                        m  = indm[q,3];
                        a1 = Ar[n,m + l + (l*l) + 1,k,sj]
                        b1 = Ai[n,m + l + (l*l) + 1,k,sj]
                        a2 = Ar[n,m1 + l1 + (l1*l1) + 1,k,sj]
                        b2 = Ai[n,m1 + l1 + (l1*l1) + 1,k,sj]
                        a3 = Ar[n,m2 + l2 + (l2*l2) + 1,k,sj]
                        b3 = Ai[n,m2 + l2 + (l2*l2) + 1,k,sj]
                        tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                         
                    end      
                    bs[n,j,k,t[n],sj] = tmp
                end
            end
        end
    end

    bs = reshape(bs, (N, J*K*S*S)) 
    if (dAr === nothing) | (dAi === nothing) 
        return bs 
    end
    
    dim, Nij, M, K = size(dAr)
    dbsa = zeros(dim,Nij,J,K);        
    Q = size(indm,1);
    bispectrumderiv(dbsa, Ar, Ai, dAr, dAi, cg, pairnumsum, indl[:], indm[:], rowm, 
        Int32(K), Int32(J), Int32(Q), Int32(M), Int32(N), Int32(Nij));

    dbs = zeros(dim,N,J,K,S,S);        
    for n = 1:Nij
        ii = ai[n]
        ij = aj[n]
        t1 = ti[n]
        t2 = tj[n]
        for k = 1:K     
            for j = 1:J
                dbs[1,ii,j,k,t1,t2] += dbsa[1,n,j,k]
                dbs[1,ij,j,k,t1,t2] -= dbsa[1,n,j,k]
                dbs[2,ii,j,k,t1,t2] += dbsa[2,n,j,k]
                dbs[2,ij,j,k,t1,t2] -= dbsa[2,n,j,k]                     
                dbs[3,ii,j,k,t1,t2] += dbsa[3,n,j,k]
                dbs[3,ij,j,k,t1,t2] -= dbsa[3,n,j,k]                     
            end
        end
    end
    dbs = reshape(dbs, (dim*N, J*K*S*S)) 

    return bs, dbs   
end

function fourbodyeatom(x, t, a, b, c, pbc, pod::PODdesc)

    rcut = pod.rcut 
    y, alist, neighlist, pairnumsum, pairnum = fullneighborlist(x, a, b, c, pbc, rcut);

    Nij = length(neighlist)        
    if Nij > 0
        dim, natom = size(x)

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

        ai = ai .+ Int32(1)
        aj = aj .+ Int32(1)    
        phi, dphi = radialsphericalharmonicbasis(pod, rij);
        A = radialsphericalharmonicbasis_sum(phi, ai, natom);
        #eatom = bispectrum(real(A), imag(A), nothing, nothing, ai, aj, Int32.(pairnumsum), pod.sh);             
        #N,J,K = size(eatom)
        #eatom = reshape(eatom, (N, J*K))
        eatom = bispectrum(real(A), imag(A), nothing, nothing, t, ai, aj, ti, tj, Int32.(pairnumsum), pod);             
        return eatom
    else
        dim, N = size(x)        
        J = size(pod.sh.indl,1);
        K = size(pod.Phi,2)
        eatom = zeros(N, J*K)
        return eatom, fatom 
    end
end

function fourbodyefatom(x, t, a, b, c, pbc, pod::PODdesc)

    rcut = pod.rcut 
    y, alist, neighlist, pairnumsum, pairnum = fullneighborlist(x, a, b, c, pbc, rcut);

    Nij = length(neighlist)        
    if Nij > 0
        dim, natom = size(x)
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

        ai = ai .+ Int32(1)
        aj = aj .+ Int32(1)                
        phi, dphi = radialsphericalharmonicbasis(pod, rij);
        A = radialsphericalharmonicbasis_sum(phi, ai, natom);
        # eatom, fatom = bispectrum(real(A), imag(A), real(dphi), imag(dphi), ai, aj, Int32.(pairnumsum), pod.sh);                    
        # dim,N,J,K = size(fatom)
        # eatom = reshape(eatom, (N, J*K))
        # fatom = reshape(fatom, (dim*N, J*K))
        eatom, fatom = bispectrum(real(A), imag(A), real(dphi), imag(dphi), t, ai, aj, ti, tj, Int32.(pairnumsum), pod);             
        return eatom, fatom 
    else
        dim, N = size(x)        
        J = size(pod.sh.indl,1);
        K = size(pod.Phi,2)
        eatom = zeros(N, J*K)
        fatom = zeros(dim*N, J*K)
        return eatom, fatom 
    end
end

function PODdescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
    if pod.nbody == 1
        return onebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod)
    elseif pod.nbody == 2
        return twobodydescriptors(x, t, a, b, c, pbc, pod)
    elseif pod.nbody == 3
        return threebodydescriptors(x, t, a, b, c, pbc, pod)
    elseif pod.nbody == 4
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
            # dim,natom,M2 = size(dd2)
            # K2 = size(pod[n].U,2)                        
            # d2 = reshape(d2,(1,M2))*pod[n].U 
            # d2 = d2[:]
            # dd2 = reshape(dd2,(dim*natom,M2))*pod[n].U             
            # dd2 = reshape(dd2, (dim,natom,K2))
            nbd2 = 1
        elseif pod[n].nbody==3
            d3, dd3 = threebodydescriptors(x, t, a, b, c, pbc, pod[n])
            # dim,natom,M3 = size(dd3)
            # K3 = size(pod[n].U,2)                        
            # d3 = reshape(d3,(1,M3))*pod[n].U 
            # d3 = d3[:]
            # dd3 = reshape(dd3,(dim*natom,M3))*pod[n].U             
            # dd3 = reshape(dd3, (dim,natom,K3))
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

    # if (nbd2==1) & (nbd3==1)
    #     dim = size(dd2,1)
    #     natom = size(dd2,2)
    #     M2 = length(d2)
    #     M3 = length(d3)

    #     d4 = zeros(M2*M3)
    #     dd4 = zeros(dim, natom, M2*M3)
    #     podhybrid23(d4, dd4, d2, d3, dd2, dd3, Int32(M2), Int32(M3), Int32(dim*natom))
    #     # for m3 = 1:M3
    #     #     for m2 = 1:M2
    #     #         d4[m2+M2*(m3-1)] = d2[m2]*d3[m3]                
    #     #         dd4[:,:,m2+M2*(m3-1)] = d2[m2]*dd3[:,:,m3] + dd2[:,:,m2]*d3[m3]
    #     #     end
    #     # end
    #     d4 = d4[:]        
    #     d = [d; d4]
    #     dd = cat(dd, dd4, dims=3)
    # end

    dim, N, M = size(dd)
    dv = zeros(6,M)

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
    B4 = 0.0
    K2 = 0
    K3 = 0
    K4 = 0
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        A2, A3, A4, N2, N3, N4 = podeatom(config, descriptors)        
        B2 = B2 .+ A2
        B3 = B3 .+ A3
        B4 = B4 .+ A4 
        K2 = K2 + N2
        K3 = K3 + N3
        K4 = K4 + N4
    end    

    if K2 > 0
        projectiontol = 1e-11
        for n = 1:length(descriptors)
            if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
                projectiontol = descriptors[n].projectiontol
                break 
            end
        end
        B2 = (1.0/K2)*B2    
        U2, s2 = eigenvalues(B2)
        n2 = numbasis(s2, projectiontol);
        U2 = U2[:,1:n2];
    end

    if K3 > 0
        projectiontol = 1e-11
        for n = 1:length(descriptors)
            if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
                projectiontol = descriptors[n].projectiontol
                break 
            end
        end
        B3 = (1.0/K3)*B3
        U3, s3 = eigenvalues(B3)    
        n3 = numbasis(s3, projectiontol);    
        U3 = U3[:,1:n3];    
    end

    if K4 > 0
        projectiontol = 1e-11
        for n = 1:length(descriptors)
            if (descriptors[n].name == "POD") & (descriptors[n].nbody==4) 
                projectiontol = descriptors[n].projectiontol
                break 
            end
        end
        B4 = (1.0/K4)*B4
        U4, s4 = eigenvalues(B4)    
        n4 = numbasis(s4, projectiontol);    
        U4 = U4[:,1:n4];    
    end

    for n = 1:length(descriptors)
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
            descriptors[n].U = U2 
        end
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
            descriptors[n].U = U3
        end
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==4) 
            descriptors[n].U = U4
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
A4 = 0.0
N2 = 0
N3 = 0 
N4 = 0
for i = 1:nconfigs
    ci = i
    x, t, a, b, c, pbc = getconfig(config, natom, ci)

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
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==4) 
            eatom = fourbodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                
            A4 = A4 .+ eatom'*eatom 
            N4 = N4 .+ size(eatom,1)
        end
    end
end

return A2, A3, A4, N2, N3, N4 

end

function PODglobaldescriptors(x, t, a, b, c, pbc, descriptors)

for n = 1:length(descriptors)
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==1) 
        eatom1, fatom1 = onebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                                        
        globd1 = sum(eatom1, dims=1)            
    end    
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
        eatom1 = twobodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                            
        globd1 = sum(eatom1, dims=1)
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
        eatom1 = threebodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                        
        globd1 = sum(eatom1, dims=1)
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==4) 
        eatom1 = fourbodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                        
        globd1 = sum(eatom1, dims=1)
    end
    if (n==1)
        globd = 1.0*globd1
    else
        globd = cat(globd, globd1, dims=2)
    end
end

globd = globd[:]

return globd
    
end
    
function podeatom(x, t, a, b, c, pbc, descriptors)

for n = 1:length(descriptors)
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==1) 
        eatom1, fatom1 = onebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                                        
        globd1 = sum(eatom1, dims=1)            
    end    
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
        eatom1 = twobodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                            
        globd1 = sum(eatom1, dims=1)
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
        eatom1 = threebodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                        
        globd1 = sum(eatom1, dims=1)
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==4) 
        eatom1 = fourbodyeatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                        
        globd1 = sum(eatom1, dims=1)
    end
    if (n==1)
        globd = 1.0*globd1
        eatom = 1.0*eatom1 
    else
        globd = cat(globd, globd1, dims=2)
        eatom = cat(eatom, eatom1, dims=2)
    end
end

globd = globd[:]
eatom = eatom 

return globd, eatom 
    
end
   

function podefatom(x, t, a, b, c, pbc, descriptors)

nbd1 = 0; nbd2 = 0; nbd3 = 0; nbd4=0;    
eatom1 = 0.0; eatom2 = 0.0; eatom3 = 0.0; eatom4 = 0.0;
fatom1 = 0.0; fatom2 = 0.0; fatom3 = 0.0; fatom4 = 0.0; fatom23 = 0.0
globd1 = 0.0; globd2 = 0.0; globd3 = 0.0; globd4 = 0.0; globd23 = 0.0
hybrid23 = false 
for n = 1:length(descriptors)
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==1) 
        eatom1, fatom1 = onebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                                        
        globd1 = sum(eatom1, dims=1)            
        nbd1 = 1
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
        eatom2, fatom2 = twobodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                            
        globd2 = sum(eatom2, dims=1)
        nbd2 = 1
        hybrid23 = descriptors[n].hybrid23
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==3) 
        eatom3, fatom3 = threebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                
        globd3 = sum(eatom3, dims=1)
        nbd3 = 1
        hybrid23 = descriptors[n].hybrid23
    end
    if (descriptors[n].name == "POD") & (descriptors[n].nbody==4) 
        eatom4, fatom4 = fourbodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                
        globd4 = sum(eatom4, dims=1)
        nbd4 = 1        
    end
end    

# if (nbd1 == 1) & (nbd2==1) & (nbd3==1)
#     globd = cat(globd1, globd2, dims=2)
#     globd = cat(globd, globd3, dims=2)
#     fatom = cat(fatom1, fatom2, dims=2)
#     fatom = cat(fatom, fatom3, dims=2)
#     eatom = cat(eatom1, eatom2, dims=2)
#     eatom = cat(eatom, eatom3, dims=2)
# elseif (nbd2==1) & (nbd3==1)
#     globd = cat(globd2, globd3, dims=2)
#     fatom = cat(fatom2, fatom3, dims=2)
#     eatom = cat(eatom2, eatom3, dims=2)
# elseif (nbd2==1) & (nbd1==1)
#     globd = cat(globd1, globd2, dims=2)
#     fatom = cat(fatom1, fatom2, dims=2)
#     eatom = cat(eatom1, eatom2, dims=2)
# elseif (nbd2==1) & (nbd3==1)
#     globd = cat(globd2, globd3, dims=2)
#     fatom = cat(fatom2, fatom3, dims=2)
#     eatom = cat(eatom2, eatom3, dims=2)
# elseif (nbd2==1) 
#     globd = 1.0*globd2         
#     fatom = 1.0*fatom2         
#     eatom = 1.0*eatom2         
# elseif (nbd3==1) 
#     globd = 1.0*globd3         
#     fatom = 1.0*fatom3         
#     eatom = 1.0*eatom3             
# end

if (nbd2==1) & (nbd1==1)
    globd = cat(globd1, globd2, dims=2)
    fatom = cat(fatom1, fatom2, dims=2)
    eatom = cat(eatom1, eatom2, dims=2)
elseif (nbd2==1) 
    globd = 1.0*globd2         
    fatom = 1.0*fatom2         
    eatom = 1.0*eatom2         
else
    error("Two-body POD descriptors are required.")
end

if (nbd3 == 1)
    globd = cat(globd, globd3, dims=2)
    fatom = cat(fatom, fatom3, dims=2)
    eatom = cat(eatom, eatom3, dims=2)
end
if (nbd4 == 1)
    globd = cat(globd, globd4, dims=2)
    fatom = cat(fatom, fatom4, dims=2)
    eatom = cat(eatom, eatom4, dims=2)
end

# if (nbd2==1) & (nbd3==1) & (hybrid23)
#     dim, N = size(x)
#     M2 = length(globd2)
#     M3 = length(globd3)

#     globd23 = zeros(1,M2*M3)
#     fatom23 = zeros(dim*N, M2*M3)    
#     podhybrid23(globd23, fatom23, globd2, globd3, fatom2, fatom3, Int32(M2), Int32(M3), Int32(dim*N))
# end           

dim, N = size(x)
n2 = 0
n3 = 0
n4 = 0
if (nbd2==1) 
    for n = 1:length(descriptors)
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==2) 
            if (length(descriptors[n].U)>1)
                globd2 = globd2*descriptors[n].U 
                fatom2 = fatom2*descriptors[n].U             
            else
                N3 = size(fatom2,1)                
                nrbf = size(descriptors[n].Phi,2)
                M2 = min(nrbf, 7)
                S = length(descriptors[n].species)
                S2 = Int32(S*(S+1)/2)
                globd2 = reshape(globd2,(nrbf, S2))
                fatom2 = reshape(fatom2,(N3, nrbf, S2))                
                globd2 = reshape(globd2[1:M2,:], (1, M2*S2))
                fatom2 = reshape(fatom2[:,1:M2,:], (N3, M2*S2))
            end
            n2 = n;
            break; 
        end
    end
end
if (nbd3==1) 
    for n = 1:length(descriptors)
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==3)  
            if (length(descriptors[n].U)>1)
                globd3 = globd3*descriptors[n].U 
                fatom3 = fatom3*descriptors[n].U          
            else
                nrbf = size(descriptors[n].Phi,2)
                nabf = descriptors[n].pdegree[3]       
                M = nrbf*(nabf+1)  
                N3 = size(fatom2,1)                
                M3 = min(nrbf, 7)*(min(nabf, 4) + 1)                
                S = length(descriptors[n].species)
                S3 = Int32(S*S*(S+1)/2)
                globd3 = reshape(globd3,(M, S3))
                fatom3 = reshape(fatom3,(N3, M, S3))                
                globd3 = reshape(globd3[1:M3,:], (1, M3*S3))
                fatom3 = reshape(fatom3[:,1:M3,:], (N3, M3*S3))
            end
            n3 = n;
            break;    
        end
    end
end
if (nbd4==1) 
    for n = 1:length(descriptors)
        if (descriptors[n].name == "POD") & (descriptors[n].nbody==4)  
            if (length(descriptors[n].U)>1)
                globd4 = globd4*descriptors[n].U 
                fatom4 = fatom4*descriptors[n].U          
            end
            n4 = n;
            break;    
        end
    end
end

M2 = descriptors[n2].b2b2
if (nbd2==1) & (M2 > 0)    
    globd22 = zeros(1,M2*M2)
    fatom22 = zeros(dim*N, M2*M2)    
    podhybrid23(globd22, fatom22, globd2, globd2, fatom2, fatom2, Int32(M2), Int32(M2), Int32(dim*N))
    globd = cat(globd, globd22, dims=2)    
    fatom = cat(fatom, fatom22, dims=2)
end

if (nbd3==1) 
    M3 = descriptors[n3].b3b3
    if (M3 > 0)    
    globd33 = zeros(1,M3*M3)
    fatom33 = zeros(dim*N, M3*M3)    
    podhybrid23(globd33, fatom33, globd3, globd3, fatom3, fatom3, Int32(M3), Int32(M3), Int32(dim*N))
    globd = cat(globd, globd33, dims=2)    
    fatom = cat(fatom, fatom33, dims=2)
    end
end

if (nbd2==1) & (nbd3==1) 
    M2 = descriptors[n2].b2b3
    M3 = descriptors[n3].b2b3
    if (M2 > 0) & (M3 > 0)
    globd23 = zeros(1,M2*M3)
    fatom23 = zeros(dim*N, M2*M3)    
    podhybrid23(globd23, fatom23, globd2, globd3, fatom2, fatom3, Int32(M2), Int32(M3), Int32(dim*N))
    globd = cat(globd, globd23, dims=2)    
    fatom = cat(fatom, fatom23, dims=2)
    end
end

if (nbd2==1) & (nbd4==1) 
    M2 = descriptors[n2].b2b4
    M4 = descriptors[n4].b2b4
    if (M2 > 0) & (M4 > 0)
    globd24 = zeros(1,M2*M4)
    fatom24 = zeros(dim*N, M2*M4)    
    podhybrid23(globd24, fatom24, globd2, globd4, fatom2, fatom4, Int32(M2), Int32(M4), Int32(dim*N))
    globd = cat(globd, globd24, dims=2)    
    fatom = cat(fatom, fatom24, dims=2)
    end
end

if (nbd3==1) & (nbd4==1) 
    M3 = descriptors[n3].b3b4
    M4 = descriptors[n4].b3b4
    if (M3 > 0) & (M4 > 0)
    globd34 = zeros(1,M3*M4)
    fatom34 = zeros(dim*N, M3*M4)    
    podhybrid23(globd34, fatom34, globd3, globd4, fatom3, fatom4, Int32(M3), Int32(M4), Int32(dim*N))
    globd = cat(globd, globd34, dims=2)    
    fatom = cat(fatom, fatom34, dims=2)
    end
end

if (nbd2==1) & (nbd3==1) & (hybrid23)
    M2 = length(globd2)
    M3 = length(globd3)
    
    globd23 = zeros(1,M2*M3)
    fatom23 = zeros(dim*N, M2*M3)    
    podhybrid23(globd23, fatom23, globd2, globd3, fatom2, fatom3, Int32(M2), Int32(M3), Int32(dim*N))

    globd = cat(globd, globd23, dims=2)    
    fatom = cat(fatom, fatom23, dims=2)
end

return globd, fatom, eatom 

end

function podlinearsystem(config, descriptors, normalizeenergy)

# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

matA = 1.0
vecb = 1.0
for i = 1:nconfigs
    ci = i
    normconst = 1.0
    if normalizeenergy==1
        normconst = 1.0/config.natom[ci]
    end    
    x, t, a, b, c, pbc, e, f = getconfig(config, natom, ci)
    globd, fatom = podefatom(x, t, a, b, c, pbc, descriptors)
    
    f = f[:];
    we2 = (config.we[1]*config.we[1])*(normconst*normconst)
    wf2 = (config.wf[1]*config.wf[1])
    if i == 1
        matA = (we2*(globd'*globd) + wf2*(fatom'*fatom))    
        vecb = ((we2*e)*(globd') + wf2*(fatom'*f))    
    else
        matA += (we2*(globd'*globd) + wf2*(fatom'*fatom))    
        vecb += ((we2*e)*(globd') + wf2*(fatom'*f))    
    end
end

return matA, vecb

end

function podleastsq(data, descriptors, Doptions)

    display("Perform linear regression fitting ...")

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
    matA = 1.0
    vecb = 1.0
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        A1, b1 = podlinearsystem(config, descriptors, normalizeenergy)        
        if i == 1
            matA = A1
            vecb = b1         
        else
            matA += A1
            vecb += b1         
        end
    end    
    
    for i = 1:size(matA,1)
        matA[i,i] = matA[i,i]*(1.0 + 1e-12);
    end

    coeff = matA\vecb 
    return coeff, matA, vecb 
end

function poderrors(config, descriptors, coeff, normalizeenergy)

    # cumalative sum of numbers of atoms 
    nconfigs = length(config.natom)
    natom = [0; cumsum(config.natom[:])];
    
    energy = zeros(nconfigs)
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

        energy[i] = e1[1]*normconst
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

    return eerr, emae, ermse, fmae, frmse, fmaeave, frmseave, nconfigs, nforces, energy
end
    
function poderroranalysis(data, descriptors, Doptions, coeff)

    display("Calculate errors ...")

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
    energies = Array{Any}(nothing, n)
    eerr = Array{Any}(nothing, n)
    ferr = Array{Any}(nothing, n)
    ferm = Array{Any}(nothing, n)
    
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        
        #ntest = config.nconfigs
        #println("Number of test configurations for " * data[i].datapath * " : $ntest")                      

        e1, e2, e3, e4, e5, e6, e7, nconfigs, nforces, energy = 
                poderrors(config, descriptors, coeff, normalizeenergy)        
        
        energies[i] = energy   
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

    return energyerrors, forceerrors, stresserrors, eerr, ferr, ferm, energies  
end

