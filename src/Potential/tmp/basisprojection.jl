function scaledbesselprojection(pdegree, dij, rin, rcut, scalefac)

    # pdegree = 6
    # rin = 1.0;
    # rcut = 5.0;
    # dij = LinRange(rin+1e-6,rcut,10000)
    # dij = collect(dij)
    # scalefac = collect(-5:0.5:5)

    epsil = 1e-6;    
    N = length(dij);
    M = length(scalefac);
    nbf = pdegree*M;
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

    rbfnew, U, s = eigendecomposition(rbf, dij)

    # A = (1.0/N)*(rbf'*rbf);
    # s = eigvals(A)
    # U = eigvecs(A)
    # ind = sortperm(s, rev=true)
    # s = s[ind]
    # U = U[:,ind]

    # # sumeig = sum(s);
    # # for i = 1:nbf
    # #     a = sum(s(1:i))/sumeig;
    # #     if a>=(1-1.0e-10)
    # #         break;
    # #     end
    # # end

    # rbfnew = rbf*U;
    # sij = dij[2:end]-dij[1:end-1];
    # # normalize the basis on [rin, rcut]
    # for i = 1:size(rbfnew,2)
    #     rbf2 = rbfnew[:,i].*rbfnew[:,i];        
    #     area = 0.5*sum(sij.*(rbf2[2:end]+rbf2[1:end-1]));    
    #     U[:,i] = U[:,i]/sqrt(area);
    # end
    # rbfnew = rbf*U; 

    return rbf, rbfnew, U, s
end

function scaledbesselprojection(pdegree, rin, rcut, scalefac)

    dij = LinRange(rin+1e-6, rcut, 2000)
    dij = collect(dij)
    #scalefac = collect(-5:0.5:5)

    return scaledbesselprojection(pdegree, dij, rin, rcut, scalefac)
end

function radialprojection(pdegree, dij, rin, rcut, scalefac)

    # pdegree = 6
    # rin = 1.0;
    # rcut = 5.0;
    # dij = LinRange(rin+1e-6,rcut,10000)
    # dij = collect(dij)
    # scalefac = collect(-5:0.5:5)

    epsil = 1e-6;    
    K = 12
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

function radialbasisprojection(pdegree, rin, rcut, scalefac)

    dij = LinRange(rin+1e-6, rcut, 2000)
    dij = collect(dij)

    rbf, rbfnew, U, s = radialprojection(pdegree, dij, rin, rcut, scalefac)
    k = numbasis(s)
    U = U[:,1:k]

    return U, s
end

function cosinprojection(pdegree, theta, theta0)
    # pdegree = 10    
    # theta = collect(LinRange(0.0, pi, 10000))     
    # theta0 = collect(LinRange(0.0, pi, 12))     
    
    N = length(theta);
    M = length(theta0);
    nbf = pdegree*M;
    abf = zeros(N, nbf);
    
    for j = 1:M    
        for i = 1:pdegree            
            abf[:,i+(j-1)*pdegree] = (cos.(theta)*cos(theta0[j]) + sin.(theta)*sin(theta0[j])).^i; 
        end
    end    

    abfnew, U, s = eigendecomposition(abf, theta)

    # A = (1.0/N)*(abf'*abf);
    # s = eigvals(A)
    # U = eigvecs(A)
    # ind = sortperm(s, rev=true)
    # s = s[ind]
    # U = U[:,ind]
    
    # abfnew = abf*U;
    # sij = theta[2:end] - theta[1:end-1];
    # # normalize the basis on [0, pi]
    # for i = 1:size(abfnew,2)
    #     abf2 = abfnew[:,i].*abfnew[:,i];        
    #     area = 0.5*sum(sij.*(abf2[2:end]+abf2[1:end-1]));    
    #     U[:,i] = U[:,i]/sqrt(area);
    # end
    # abfnew = abf*U; 

    return abf, abfnew, U, s
end

function cosinbasisprojection(pdegree, theta0)

    theta = collect(LinRange(0.0, pi, 2000))     
    rbf, rbfnew, U, s = cosinprojection(pdegree, theta, theta0)
    k = numbasis(s)
    U = U[:,1:k]

    return U, s
end

function fourierprojection(pdegree, theta, theta0)
    # pdegree = 10    
    # theta = collect(LinRange(0.0, pi, 10000))     
    # theta0 = collect(LinRange(0.0, pi, 12))     
    
    N = length(theta);
    M = length(theta0);
    nbf = pdegree*M;
    abf = zeros(N, nbf);
    
    for j = 1:M    
        for i = 1:pdegree            
            abf[:,i+(j-1)*pdegree] = cos.(i*theta .- theta0[j]) 
        end
    end    

    abfnew, U, s = eigendecomposition(abf, theta)

    return abf, abfnew, U, s
end

function fourierbasisprojection(pdegree, theta0)

    theta = collect(LinRange(0.0, pi, 2000))     
    rbf, rbfnew, U, s = fourierprojection(pdegree, theta, theta0)
    k = numbasis(s)
    U = U[:,1:k]

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

function numbasis(s, tol=1e-6)
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












