

@kwdef mutable struct PODdesc 
    name::String = "POD"
    nbody::Int64 = 2
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
    atomtypes::Vector{Int64} = [1, 1, 1, 1], # a list of atom types 
    pdegree::Vector{Int64} = [5,12,5], # radial and angular polynomial degrees 
    nbasis::Vector{Int64} = [5,5], # number of radial and angular descriptors 
    rcut::Float64 = 5.0,    # cut-off radius 
    rin::Float64 = 1.0,     # inner cut-off radius 
    gamma0::Vector{Float64} = [-2, 0.0, 2], # radial parameters
    theta0::Vector{Float64} = [0.0],    # angular parameters
    Phi::Matrix{Float64} = zeros(1,1), # radial projection matrix 
    Psi::Matrix{Float64} = zeros(1,1), # angular projection matrix 
    Lambda::Vector{Float64} = [0.0],   # radial eigenvalues
    Upsilon::Vector{Float64} = [0.0],  # angular eigenvalues
    )
        
    pod = initPOD(PODdesc(name, nbody, atomtypes, pdegree, nbasis, rcut, rin, gamma0, theta0, Phi, Psi, Lambda, Upsilon))

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
    # k = numbasis(s, tol)
    # U = U[:,1:k]

    return U, s
end

# function fourierprojection(pdegree, theta, theta0)
    
#     N = length(theta);
#     M = length(theta0);
#     nbf = pdegree*M;
#     abf = zeros(N, nbf);
    
#     for j = 1:M    
#         for i = 1:pdegree            
#             abf[:,i+(j-1)*pdegree] = cos.(i*theta .+ theta0[j]) 
#         end
#     end    

#     abfnew, U, s = eigendecomposition(abf, theta)

#     return abf, abfnew, U, s
# end

# function fourierbasisprojection(pdegree, theta0, tol)

#     theta = collect(LinRange(0.0, pi, 2000))     
#     rbf, rbfnew, U, s = fourierprojection(pdegree, theta, theta0)
#     # k = numbasis(s, tol)
#     # U = U[:,1:k]

#     return U, s 
# end

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

function radialbasis(rij, rin, rcut, scalefac, pdegree, K)

    epsil = 1e-6;    
    dim, N = size(rij);
    if dim==3
        dij = sqrt.(rij[1,:].*rij[1,:] .+ rij[2,:].*rij[2,:] .+ rij[3,:].*rij[3,:])    
    else
        dij = sqrt.(rij[1,:].*rij[1,:] .+ rij[2,:].*rij[2,:])
    end
    dr = rij./reshape(dij,(1,N));
    
    M = length(scalefac);
    nbf = pdegree*M + K
    rbf = zeros(N, nbf);
    drbf = zeros(dim, N, nbf);

    r = dij .- rin
    rmax = rcut - rin
    #fcut = cutoff(dij, rin, rcut, epsil) 
    fcut,dfcut = cutoffderiv(dij, rin, rcut, epsil) 

    for j = 1:M
        alpha = scalefac[j];    
        if alpha == 0
            alpha = 1e-3;
        end
        x =  (1.0 .- exp.(-alpha*r/rmax))/(1.0-exp(-alpha));
        dx = (alpha/rmax)*exp.(-(alpha*r/rmax))/(1.0 - exp(-alpha));        

        for i = 1:pdegree    
            a = i*pi;
            rbf[:,i+(j-1)*pdegree] = (sqrt(2.0/(rmax))/i)*fcut.*sin.(a*x)./r;
            drbfdr = (sqrt(2/rmax)/i)*(dfcut.*sin.(a*x)./r - fcut.*sin.(a*x)./(r.^2) + a*cos.(a*x).*fcut.*dx./r);
            drbf[:,:,i+(j-1)*pdegree] = reshape(drbfdr, (1, N)).*dr;
            # b = (sqrt(2/rcut)/i)*(sin.(a*x)./(dij.^2) .- (a*cos.(a*x).*dx)./dij);
            # drbf[:,:,i] = -dr.*reshape(b,(1,N))        
        end
    end

    # add radial basis
    for i = 1:K
        n = pdegree*M+i;
        rbf[:,n] = fcut./(dij.^i);
        drbfdr = dfcut./(dij.^i) - i*fcut./(dij.^(i+1));  
        drbf[:,:,n] = reshape(drbfdr, (1, N)).*dr;
    end
    
    return rbf, drbf 
end

function radialbasis2(rij, rin, rcut, scalefac, pdegree, K)

    epsil = 1e-6;    
    dim, N = size(rij);
    M = length(scalefac);
    nbf = pdegree*M + K
    rbf = zeros(N, nbf);
    drbf = zeros(dim, N, nbf);
    
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
                drbf[1,n,i+(j-1)*pdegree] = drbfdr*dr1;
                drbf[2,n,i+(j-1)*pdegree] = drbfdr*dr2;
                drbf[3,n,i+(j-1)*pdegree] = drbfdr*dr3;
            end
        end

        # add radial basis
        for i = 1:K
            p = pdegree*M+i;
            rbf[n,p] = fcut/(dij^i);
            drbfdr = dfcut/(dij^i) - i*fcut/(dij^(i+1));  
            #drbf[:,n,p] = drbfdr*dr;
            drbf[1,n,p] = drbfdr*dr1;
            drbf[2,n,p] = drbfdr*dr2;
            drbf[3,n,p] = drbfdr*dr3;
        end
    end
    
    return rbf, drbf 
end

function onebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, N = size(x)
    typei = pod.atomtypes[1];    
    ilist = findatomtype(Array(1:N), t, typei);      
    xi = x[:,ilist]
    ai = ilist
    ti = t[ilist]    
    ei = ones(length(ti))
    fi = zeros(size(xi))        
    ea, fa = tallysingle(ei, fi, ai);   
    d = [sum(ea)]         
    dd = fa
    va = tallysinglevirial(fi, xi, ai, -1.0)
    dv = sum(va, dims=2)
    dd = reshape(dd, (dim,size(fa,2),1))        
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
    M = size(pod.Phi,2)
    for n = 1:dim 
        fij[n,:,1:M] = fij[n,:,:]*pod.Phi 
    end

    d = []
    dd = []
    dv = []
    for m = 1:M
        ea, fa = tallypair(eij[:,m], fij[:,:,m], ai, aj);            
        d = [d; sum(ea)];
        if length(dd)==0
            dd = fa;
        else
            dd = cat(dd, fa, dims=3)
        end    
        va = tallypairvirial(fij[:,:,m], rij, ai, aj, -1.0/2)
        if length(dv)==0
            dv = sum(va, dims=2)
        else
            dv = [dv sum(va, dims=2)]
        end                
    end
    
    return d, dd, dv 
end

function fourierbasis(xij, xik, nabf)

    dim, N = size(xij);

    xij1 = xij[1,:]
    xij2 = xij[2,:]
    xij3 = xij[3,:]
    xik1 = xik[1,:]
    xik2 = xik[2,:]
    xik3 = xik[3,:]

    xdot  = xij1.*xik1 + xij2.*xik2 + xij3.*xik3;
    rijsq = xij1.*xij1 + xij2.*xij2 + xij3.*xij3;    
    riksq = xik1.*xik1 + xik2.*xik2 + xik3.*xik3;    
    rij = sqrt.(rijsq); 
    rik = sqrt.(riksq); 

    costhe = xdot./(rij.*rik);    
    ind = findall(abs.(costhe) .>= 1.0)
    costhe[ind] = 1.0*sign.(costhe[ind])
    xdot = costhe.*(rij.*rik)

    sinthe = sqrt.(1 .- costhe.*costhe);
    sinthe = max.(sinthe, 1e-12)    
    theta = acos.(costhe)            
    dtheta = -1.0./sinthe 

    tm1 = (rijsq.^(3/2)).*(riksq.^(1/2));
    tm2 = (rijsq.^(1/2)).*(riksq.^(3/2));
    dct1 = (xik1.*rijsq - xij1.*xdot)./tm1; 
    dct2 = (xik2.*rijsq - xij2.*xdot)./tm1;
    dct3 = (xik3.*rijsq - xij3.*xdot)./tm1;
    dct4 = (xij1.*riksq - xik1.*xdot)./tm2;
    dct5 = (xij2.*riksq - xik2.*xdot)./tm2;
    dct6 = (xij3.*riksq - xik3.*xdot)./tm2;

    abf  = zeros(N, 2*nabf)
    dabf  = zeros(6, N, 2*nabf)    
    for p = 1:nabf
        abf[:,p] = cos.(p*theta);                
        tm = -p*sin.(p*theta).*dtheta
        dabf[1,:,p] = tm.*dct1;
        dabf[2,:,p] = tm.*dct2;
        dabf[3,:,p] = tm.*dct3;
        dabf[4,:,p] = tm.*dct4;
        dabf[5,:,p] = tm.*dct5;
        dabf[6,:,p] = tm.*dct6;
        q = nabf+p
        abf[:,q] = sin.(p*theta);                
        tm = p*cos.(p*theta).*dtheta
        dabf[1,:,q] = tm.*dct1;
        dabf[2,:,q] = tm.*dct2;
        dabf[3,:,q] = tm.*dct3;
        dabf[4,:,q] = tm.*dct4;
        dabf[5,:,q] = tm.*dct5;
        dabf[6,:,q] = tm.*dct6;
    end
        
    return abf, dabf 
end

function cosinbasis(xij, xik, nabf, theta0)

    dim, N = size(xij);
    n0 = length(theta0)
    abf  = zeros(N, n0*nabf)
    dabf  = zeros(6, N, n0*nabf)        

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

    return abf, dabf 
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
    nn = size(e2ij,1)
    f2ij = reshape(reshape(f2ij,(3*nn,nm))*pod.Phi, (3, nn, nrbf)) 
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

    #@time begin
    nabf = size(uijk,2)
    M = nrbf*(nabf+1)  
    N = size(uij,1)    
    fij = zeros(3)
    fik = zeros(3)    
    natom = size(x,2)
    eatom = zeros(natom,M)
    fatom = zeros(3,natom,M)
    vatom = zeros(6,natom,M)    
    for n = 1:N
        K = 0
        for m = 1:(nrbf)
            K = K + 1
            eijk = uij[n,m]*uik[n,m];
            fij[1] = wij[1,n,m]*uik[n,m];
            fij[2] = wij[2,n,m]*uik[n,m];
            fij[3] = wij[3,n,m]*uik[n,m];
            fik[1] = uij[n,m]*wik[1,n,m];
            fik[2] = uij[n,m]*wik[2,n,m];
            fik[3] = uij[n,m]*wik[3,n,m];            

            v0 = xij[1,n]*fij[1] + xik[1,n]*fik[1];
            v1 = xij[2,n]*fij[2] + xik[2,n]*fik[2];        
            v2 = xij[3,n]*fij[3] + xik[3,n]*fik[3];
            v3 = xij[2,n]*fij[3] + xik[2,n]*fik[3];
            v4 = xij[1,n]*fij[3] + xik[1,n]*fik[3];            
            v5 = xij[1,n]*fij[2] + xik[1,n]*fik[2];                                

            i1 = ai[n]
            eatom[i1,K] += eijk
            fatom[1,i1,K] += fij[1] + fik[1]
            fatom[2,i1,K] += fij[2] + fik[2]
            fatom[3,i1,K] += fij[3] + fik[3]
            vatom[1,i1,K] += v0; 
            vatom[2,i1,K] += v1;
            vatom[3,i1,K] += v2; 
            vatom[4,i1,K] += v3;
            vatom[5,i1,K] += v4; 
            vatom[6,i1,K] += v5;        
            
            i1 = aj[n]
            fatom[1,i1,K] -= fij[1]
            fatom[2,i1,K] -= fij[2]
            fatom[3,i1,K] -= fij[3]
            vatom[1,i1,K] += v0; 
            vatom[2,i1,K] += v1;
            vatom[3,i1,K] += v2; 
            vatom[4,i1,K] += v3;
            vatom[5,i1,K] += v4; 
            vatom[6,i1,K] += v5;        

            i1 = ak[n]
            fatom[1,i1,K] -= fik[1]   
            fatom[2,i1,K] -= fik[2]   
            fatom[3,i1,K] -= fik[3]   
            vatom[1,i1,K] += v0; 
            vatom[2,i1,K] += v1;
            vatom[3,i1,K] += v2; 
            vatom[4,i1,K] += v3;
            vatom[5,i1,K] += v4; 
            vatom[6,i1,K] += v5;        

            for p = 1:(nabf)                
                K = K + 1                
                # eijk = uij[n,m]*uik[n,m]*uijk[n,p];
                # fij[1] = wij[1,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[1,n,p];
                # fij[2] = wij[2,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[2,n,p];
                # fij[3] = wij[3,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[3,n,p];
                # fik[1] = uij[n,m]*wik[1,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[4,n,p];
                # fik[2] = uij[n,m]*wik[2,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[5,n,p];
                # fik[3] = uij[n,m]*wik[3,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[6,n,p];
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
                
                v0 = xij[1,n]*fij[1] + xik[1,n]*fik[1];
                v1 = xij[2,n]*fij[2] + xik[2,n]*fik[2];        
                v2 = xij[3,n]*fij[3] + xik[3,n]*fik[3];
                v3 = xij[2,n]*fij[3] + xik[2,n]*fik[3];
                v4 = xij[1,n]*fij[3] + xik[1,n]*fik[3];            
                v5 = xij[1,n]*fij[2] + xik[1,n]*fik[2];                                
    
                i1 = ai[n]
                eatom[i1,K] += eijk
                fatom[1,i1,K] += fij[1] + fik[1]
                fatom[2,i1,K] += fij[2] + fik[2]
                fatom[3,i1,K] += fij[3] + fik[3]
                vatom[1,i1,K] += v0; 
                vatom[2,i1,K] += v1;
                vatom[3,i1,K] += v2; 
                vatom[4,i1,K] += v3;
                vatom[5,i1,K] += v4; 
                vatom[6,i1,K] += v5;        
                
                i1 = aj[n]
                fatom[1,i1,K] -= fij[1]
                fatom[2,i1,K] -= fij[2]
                fatom[3,i1,K] -= fij[3]
                vatom[1,i1,K] += v0; 
                vatom[2,i1,K] += v1;
                vatom[3,i1,K] += v2; 
                vatom[4,i1,K] += v3;
                vatom[5,i1,K] += v4; 
                vatom[6,i1,K] += v5;        
    
                i1 = ak[n]
                fatom[1,i1,K] -= fik[1]   
                fatom[2,i1,K] -= fik[2]   
                fatom[3,i1,K] -= fik[3]   
                vatom[1,i1,K] += v0; 
                vatom[2,i1,K] += v1;
                vatom[3,i1,K] += v2; 
                vatom[4,i1,K] += v3;
                vatom[5,i1,K] += v4; 
                vatom[6,i1,K] += v5;                            
            end
        end
    end
    #end

    # ai = Int32.(ai.-1)
    # aj = Int32.(aj.-1)
    # ak = Int32.(ak.-1)    
    # podtally3(eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, 
    #         ai, aj, ak, Int32(nrbf), Int32(nabf), Int32(natom), Int32(N));
        
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
    # podtally3(eatom, fatom, vatom, xij, xik, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
    #     nrbf, nabf, natom, N);

    d = sum(eatom, dims=1)
    d = d[:]
    dd = fatom 
    dv = (-1.0/3.0)*reshape(sum(vatom, dims=2), (6,M))

    return d, dd, dv 
end

# function angulardescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
#     dim, N = size(x)
#     rcut = pod.rcut 
#     y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);

#     typei = pod.atomtypes[1];
#     typej = pod.atomtypes[2];            
#     typek = pod.atomtypes[3];            
#     ilist = findatomtype(Array(1:N), t, typei);      
#     pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej; typek], neighlist, neighnum, rcut*rcut);                                                
#     # ilist = Array(1:N);                     
#     # pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);     

#     #@time begin
#     tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
#     xij, xik, ai, aj, ak = neightriplets(y, [], tripletlist, tripletnum, t, ilist, alist);
#     #end

#     #@time begin
#     rij,~ = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);                
#     e2ij, f2ij = radialbasis(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
#     e2ij = e2ij*pod.Phi 
#     nm, nrbf = size(pod.Phi)    
#     nn = size(e2ij,1)
#     f2ij = reshape(reshape(f2ij,(3*nn,nm))*pod.Phi, (3, nn, nrbf)) 
#     uij, uik, wij, wik = makejk(e2ij, f2ij, pairnum, tripletnum, ilist)
#     #end

#     # @time begin
#     # uij, wij = radialbasis(xij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
#     # uik, wik = radialbasis(xik, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
#     # end
#     #@time begin
#     #uijk, wijk = fourierbasis(xij, xik, pod.pdegree[3])
#     uijk, wijk = cosinbasis(xij, xik, pod.pdegree[3], pod.theta0)
#     #end

#     #@time begin
#     # uij = uij*pod.Phi 
#     # uik = uik*pod.Phi     
#     # nm, nrbf = size(pod.Phi)
#     # nn = size(wij,2)
#     # wij = reshape(reshape(wij,(3*nn,nm))*pod.Phi, (3, nn, nrbf)) 
#     # wik = reshape(reshape(wik,(3*nn,nm))*pod.Phi, (3, nn, nrbf)) 
#     # for n = 1:dim 
#     #     wij[n,:,1:nrbf] = wij[n,:,:]*pod.Phi 
#     #     wik[n,:,1:nrbf] = wik[n,:,:]*pod.Phi 
#     # end
#     #end
#     # uijk = uijk*pod.Psi 
#     # nabf = size(pod.Psi,2)
#     # for n = 1:(2*dim) 
#     #     wijk[n,:,1:nabf] = wijk[n,:,:]*pod.Psi 
#     # end
    
#     # display(size(uij))
#     # display(size(uij2))
#     # display(maximum(abs.(uij2[:]-uij[:])))
#     # display(maximum(abs.(wij2[:]-wij[:])))
#     # display(maximum(abs.(uik2[:]-uik[:])))
#     # display(maximum(abs.(wik2[:]-wik[:])))

#     # nabf = size(uijk,2)

#     # M = nrbf*(nabf+1)  
#     # N = size(uij,1)
#     # eijk = zeros(N,M)
#     # fij = zeros(3,N,M)
#     # fik = zeros(3,N,M)    
#     # K = 0
#     # # @time begin        
#     # for m = 1:(nrbf)
#     #     K = K + 1
#     #     eijk[:,K] = uij[:,m].*uik[:,m];
#     #     fij[1,:,K] = wij[1,:,m].*uik[:,m];
#     #     fij[2,:,K] = wij[2,:,m].*uik[:,m];
#     #     fij[3,:,K] = wij[3,:,m].*uik[:,m];
#     #     fik[1,:,K] = uij[:,m].*wik[1,:,m];
#     #     fik[2,:,K] = uij[:,m].*wik[2,:,m];
#     #     fik[3,:,K] = uij[:,m].*wik[3,:,m];
#     #     for p = 1:(nabf)                
#     #         K = K + 1                
#     #         eijk[:,K] = uij[:,m].*uik[:,m].*uijk[:,p];
#     #         fij[1,:,K] = wij[1,:,m].*uik[:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[1,:,p];
#     #         fij[2,:,K] = wij[2,:,m].*uik[:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[2,:,p];
#     #         fij[3,:,K] = wij[3,:,m].*uik[:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[3,:,p];
#     #         fik[1,:,K] = uij[:,m].*wik[1,:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[4,:,p];
#     #         fik[2,:,K] = uij[:,m].*wik[2,:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[5,:,p];
#     #         fik[3,:,K] = uij[:,m].*wik[3,:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[6,:,p];
#     #     end
#     # end
#     # end
#     #@time begin        
#     # for n = 1:N
#     #     K = 0
#     #     for m = 1:(nrbf)
#     #         K = K + 1
#     #         eijk[n,K] = uij[n,m]*uik[n,m];
#     #         fij[1,n,K] = wij[1,n,m]*uik[n,m];
#     #         fij[2,n,K] = wij[2,n,m]*uik[n,m];
#     #         fij[3,n,K] = wij[3,n,m]*uik[n,m];
#     #         fik[1,n,K] = uij[n,m]*wik[1,n,m];
#     #         fik[2,n,K] = uij[n,m]*wik[2,n,m];
#     #         fik[3,n,K] = uij[n,m]*wik[3,n,m];
#     #         for p = 1:(nabf)                
#     #             K = K + 1                
#     #             eijk[n,K] = uij[n,m]*uik[n,m]*uijk[n,p];
#     #             fij[1,n,K] = wij[1,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[1,n,p];
#     #             fij[2,n,K] = wij[2,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[2,n,p];
#     #             fij[3,n,K] = wij[3,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[3,n,p];
#     #             fik[1,n,K] = uij[n,m]*wik[1,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[4,n,p];
#     #             fik[2,n,K] = uij[n,m]*wik[2,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[5,n,p];
#     #             fik[3,n,K] = uij[n,m]*wik[3,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[6,n,p];
#     #         end
#     #     end
#     # end
#     #end    
#     #error("here")

#     # d = []
#     # dd = []
#     # dv = []
#     # #@time begin
#     # for m = 1:M
#     #     #@time begin  
#     #     ea, fa = tallytriplet(eijk[:,m], fij[:,:,m], fik[:,:,m], ai, aj, ak);                                             
#     #     #end
#     #     d = [d; sum(ea)];
#     #     if length(dd)==0
#     #         dd = fa;
#     #     else
#     #         dd = cat(dd, fa, dims=3)
#     #     end    
#     #     #@time begin  
#     #     va = tallytripletvirial(fij[:,:,m], fik[:,:,m], xij, xik, ai, aj, ak, -1.0/3)
#     #     #end
#     #     if length(dv)==0
#     #         dv = sum(va, dims=2)
#     #     else
#     #         dv = [dv sum(va, dims=2)]
#     #     end                
#     # end              
#     # # end

#     # natom = size(x,2)
#     # eatom = zeros(natom,M)
#     # fatom = zeros(3,natom,M)
#     # vatom = zeros(6,natom,M)
#     # #@time begin
#     # for m = 1:M
#     #     for i = 1:N
#     #         eatom[ai[i],m] += eijk[i,m]

#     #         v0 = xij[1,i]*fij[1,i,m] + xik[1,i]*fik[1,i,m];
#     #         v1 = xij[2,i]*fij[2,i,m] + xik[2,i]*fik[2,i,m];        
#     #         v2 = xij[3,i]*fij[3,i,m] + xik[3,i]*fik[3,i,m];
#     #         v3 = xij[2,i]*fij[3,i,m] + xik[2,i]*fik[3,i,m];
#     #         v4 = xij[1,i]*fij[3,i,m] + xik[1,i]*fik[3,i,m];            
#     #         v5 = xij[1,i]*fij[2,i,m] + xik[1,i]*fik[2,i,m];                                

#     #         i1 = ai[i]
#     #         fatom[1,i1,m] += fij[1,i,m] + fik[1,i,m]
#     #         fatom[2,i1,m] += fij[2,i,m] + fik[2,i,m]
#     #         fatom[3,i1,m] += fij[3,i,m] + fik[3,i,m]
#     #         vatom[1,i1,m] += v0; 
#     #         vatom[2,i1,m] += v1;
#     #         vatom[3,i1,m] += v2; 
#     #         vatom[4,i1,m] += v3;
#     #         vatom[5,i1,m] += v4; 
#     #         vatom[6,i1,m] += v5;        
            
#     #         i1 = aj[i]
#     #         fatom[1,i1,m] -= fij[1,i,m]
#     #         fatom[2,i1,m] -= fij[2,i,m]
#     #         fatom[3,i1,m] -= fij[3,i,m]
#     #         vatom[1,i1,m] += v0; 
#     #         vatom[2,i1,m] += v1;
#     #         vatom[3,i1,m] += v2; 
#     #         vatom[4,i1,m] += v3;
#     #         vatom[5,i1,m] += v4; 
#     #         vatom[6,i1,m] += v5;        

#     #         i1 = ak[i]
#     #         fatom[1,i1,m] -= fik[1,i,m]   
#     #         fatom[2,i1,m] -= fik[2,i,m]   
#     #         fatom[3,i1,m] -= fik[3,i,m]   
#     #         vatom[1,i1,m] += v0; 
#     #         vatom[2,i1,m] += v1;
#     #         vatom[3,i1,m] += v2; 
#     #         vatom[4,i1,m] += v3;
#     #         vatom[5,i1,m] += v4; 
#     #         vatom[6,i1,m] += v5;        

#     #         # for d = 1:3
#     #         #     fatom[d,ai[i],m] += fij[d,i,m] + fik[d,i,m]
#     #         #     fatom[d,aj[i],m] -= fij[d,i,m]
#     #         #     fatom[d,ak[i],m] -= fik[d,i,m]                
#     #         # end
#     #     end
#     # end
#     # #end

#     nabf = size(uijk,2)
#     M = nrbf*(nabf+1)  
#     N = size(uij,1)    
#     fij = zeros(3)
#     fik = zeros(3)    
#     natom = size(x,2)
#     eatom = zeros(natom,M)
#     fatom = zeros(3,natom,M)
#     vatom = zeros(6,natom,M)
#     for n = 1:N
#         K = 0
#         for m = 1:(nrbf)
#             K = K + 1
#             eijk = uij[n,m]*uik[n,m];
#             fij[1] = wij[1,n,m]*uik[n,m];
#             fij[2] = wij[2,n,m]*uik[n,m];
#             fij[3] = wij[3,n,m]*uik[n,m];
#             fik[1] = uij[n,m]*wik[1,n,m];
#             fik[2] = uij[n,m]*wik[2,n,m];
#             fik[3] = uij[n,m]*wik[3,n,m];            

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

#             for p = 1:(nabf)                
#                 K = K + 1                
#                 eijk = uij[n,m]*uik[n,m]*uijk[n,p];
#                 fij[1] = wij[1,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[1,n,p];
#                 fij[2] = wij[2,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[2,n,p];
#                 fij[3] = wij[3,n,m]*uik[n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[3,n,p];
#                 fik[1] = uij[n,m]*wik[1,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[4,n,p];
#                 fik[2] = uij[n,m]*wik[2,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[5,n,p];
#                 fik[3] = uij[n,m]*wik[3,n,m]*uijk[n,p] + uij[n,m]*uik[n,m]*wijk[6,n,p];
           
#                 v0 = xij[1,n]*fij[1] + xik[1,n]*fik[1];
#                 v1 = xij[2,n]*fij[2] + xik[2,n]*fik[2];        
#                 v2 = xij[3,n]*fij[3] + xik[3,n]*fik[3];
#                 v3 = xij[2,n]*fij[3] + xik[2,n]*fik[3];
#                 v4 = xij[1,n]*fij[3] + xik[1,n]*fik[3];            
#                 v5 = xij[1,n]*fij[2] + xik[1,n]*fik[2];                                
    
#                 i1 = ai[n]
#                 eatom[i1,K] += eijk
#                 fatom[1,i1,K] += fij[1] + fik[1]
#                 fatom[2,i1,K] += fij[2] + fik[2]
#                 fatom[3,i1,K] += fij[3] + fik[3]
#                 vatom[1,i1,K] += v0; 
#                 vatom[2,i1,K] += v1;
#                 vatom[3,i1,K] += v2; 
#                 vatom[4,i1,K] += v3;
#                 vatom[5,i1,K] += v4; 
#                 vatom[6,i1,K] += v5;        
                
#                 i1 = aj[n]
#                 fatom[1,i1,K] -= fij[1]
#                 fatom[2,i1,K] -= fij[2]
#                 fatom[3,i1,K] -= fij[3]
#                 vatom[1,i1,K] += v0; 
#                 vatom[2,i1,K] += v1;
#                 vatom[3,i1,K] += v2; 
#                 vatom[4,i1,K] += v3;
#                 vatom[5,i1,K] += v4; 
#                 vatom[6,i1,K] += v5;        
    
#                 i1 = ak[n]
#                 fatom[1,i1,K] -= fik[1]   
#                 fatom[2,i1,K] -= fik[2]   
#                 fatom[3,i1,K] -= fik[3]   
#                 vatom[1,i1,K] += v0; 
#                 vatom[2,i1,K] += v1;
#                 vatom[3,i1,K] += v2; 
#                 vatom[4,i1,K] += v3;
#                 vatom[5,i1,K] += v4; 
#                 vatom[6,i1,K] += v5;                            
#             end
#         end
#     end

#     d = sum(eatom, dims=1)
#     d = d[:]
#     dd = fatom 
#     dv = (-1.0/3.0)*reshape(sum(vatom, dims=2), (6,M))

#     # display(size(dv1))
#     # display(size(dv))
#     # display(maximum(abs.(d1[:]-d[:])))
#     # display(maximum(abs.(dd1[:]-dd[:])))
#     # display(maximum(abs.(dv1[:]-dv[:])))
#     #error("here")      
#     return d, dd, dv 
# end

function angulardescriptors2(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
        
    dim, N = size(x)
    rcut = pod.rcut 
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);

    typei = pod.atomtypes[1];
    typej = pod.atomtypes[2];            
    typek = pod.atomtypes[3];            
    ilist = findatomtype(Array(1:N), t, typei);      
    pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej; typek], neighlist, neighnum, rcut*rcut);                                                
    # ilist = Array(1:N);                     
    # pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);     

    tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
    xij, xik, ai, aj, ak = neightriplets(y, [], tripletlist, tripletnum, t, ilist, alist);

    # rij,~ = neighpairs(y, [], pairlist, pairnum, t, ilist, alist);                
    # e2ij, f2ij = radialbasis(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    # e2ij = e2ij*pod.Phi 
    # M = size(pod.Phi,2)
    # for n = 1:dim 
    #     f2ij[n,:,1:M] = f2ij[n,:,:]*pod.Phi 
    # end    
    
    uij, wij = radialbasis(xij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    uik, wik = radialbasis(xik, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    #uijk, wijk = fourierbasis(xij, xik, pod.pdegree[3])
    uijk, wijk = cosinbasis(xij, xik, pod.pdegree[3], pod.theta0)

    uij = uij*pod.Phi 
    uik = uik*pod.Phi 
    nrbf = size(pod.Phi,2)
    for n = 1:dim 
        wij[n,:,1:nrbf] = wij[n,:,:]*pod.Phi 
        wik[n,:,1:nrbf] = wik[n,:,:]*pod.Phi 
    end
    # uijk = uijk*pod.Psi 
    # nabf = size(pod.Psi,2)
    # for n = 1:(2*dim) 
    #     wijk[n,:,1:nabf] = wijk[n,:,:]*pod.Psi 
    # end
    
    nabf = size(uijk,2)

    M = nrbf*(nabf+1)  
    N = size(uij,1)
    eijk = zeros(N,M)
    fij = zeros(3,N,M)
    fik = zeros(3,N,M)    
    K = 0
    for m = 1:(nrbf)
        K = K + 1
        eijk[:,K] = uij[:,m].*uik[:,m];
        fij[1,:,K] = wij[1,:,m].*uik[:,m];
        fij[2,:,K] = wij[2,:,m].*uik[:,m];
        fij[3,:,K] = wij[3,:,m].*uik[:,m];
        fik[1,:,K] = uij[:,m].*wik[1,:,m];
        fik[2,:,K] = uij[:,m].*wik[2,:,m];
        fik[3,:,K] = uij[:,m].*wik[3,:,m];
        for p = 1:(nabf)                
            K = K + 1                
            eijk[:,K] = uij[:,m].*uik[:,m].*uijk[:,p];
            fij[1,:,K] = wij[1,:,m].*uik[:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[1,:,p];
            fij[2,:,K] = wij[2,:,m].*uik[:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[2,:,p];
            fij[3,:,K] = wij[3,:,m].*uik[:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[3,:,p];
            fik[1,:,K] = uij[:,m].*wik[1,:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[4,:,p];
            fik[2,:,K] = uij[:,m].*wik[2,:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[5,:,p];
            fik[3,:,K] = uij[:,m].*wik[3,:,m].*uijk[:,p] + uij[:,m].*uik[:,m].*wijk[6,:,p];
        end
    end

    d = []
    dd = []
    dv = []    
    for m = 1:size(eijk,2)        
        ea, fa = tallytriplet(eijk[:,m], fij[:,:,m], fik[:,:,m], ai, aj, ak);                                             
        d = [d; sum(ea)];
        if length(dd)==0
            dd = fa;
        else
            dd = cat(dd, fa, dims=3)
        end    
        va = tallytripletvirial(fij[:,:,m], fik[:,:,m], xij, xik, ai, aj, ak, -1.0/3)
        if length(dv)==0
            dv = sum(va, dims=2)
        else
            dv = [dv sum(va, dims=2)]
        end                
    end                
    return d, dd, dv 
end

function PODdescriptors(x, t, a, b, c, pbc, rcutmax, pod::PODdesc)
    if pod.nbody == 1
        return onebodydescriptors(x, t, a, b, c, pbc, rcutmax, pod)
    elseif pod.nbody == 2
        # @time begin
        # radialdescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        # end
        return radialdescriptors(x, t, a, b, c, pbc, rcutmax, pod)
    elseif pod.nbody == 3
        # @time begin
        # angulardescriptors(x, t, a, b, c, pbc, rcutmax, pod)
        # end
        return angulardescriptors(x, t, a, b, c, pbc, rcutmax, pod)
    else
        error("nbody is invalid")        
    end    
end

