function podradialbasis(rij, rin, rcut, scalefac, pdegree, K)

    epsil = 1e-6;    
    dim, N = size(rij);
    M = length(scalefac);
    nbf = pdegree*M + K
    rbf = zeros(N, nbf);
    drbf = zeros(dim*N, nbf);
    
    rmax = rcut - rin
    for n = 1:N
        dij = sqrt(rij[1,n]*rij[1,n] + rij[2,n]*rij[2,n] + rij[3,n]*rij[3,n])    
        dr1 = rij[1,n]/dij;    
        dr2 = rij[2,n]/dij;    
        dr3 = rij[3,n]/dij;    
        r = dij - rin        
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
            drbf[1+3*(n-1),p] = drbfdr*dr1;
            drbf[2+3*(n-1),p] = drbfdr*dr2;
            drbf[3+3*(n-1),p] = drbfdr*dr3;
        end
    end
    
    #drbf = reshape(drbf, (dim,N,nbf));
    return rbf, drbf 
end

function podradialbasis(pod, rij::Matrix{Float64})
    eij, fij = podradialbasis(rij, pod.rin, pod.rcut, pod.gamma0, pod.pdegree[1], pod.pdegree[2])
    eij = eij*pod.Phi 
    fij = fij*pod.Phi 
    return eij, fij 
end

function radialsphericalharmonicbasis(pod, sh, rij::Matrix{Float64})

    g, ~ = podradialbasis(pod, rij)
    Ylm = cYlm!(sh, rij)

    N, K = size(g)
    M = size(Ylm,2)
    phi = zeros(Complex{Float64},N,M,K);
    for k = 1:K
        for m = 1:M    
            phi[:,m,k] = Ylm[:,m].*g[:,k]
        end
    end

    return phi
end

function radialsphericalharmonicbasis_de(pod, sh, rij)

    g, dg = podradialbasis(pod, rij)    
    Ylm, dYlm = cYlm_de!(sh, rij)

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
end

function powerspectrum1(Ar, Ai, dAr, dAi, ai, aj)

    N, M, K = size(Ar)
    L = Int32(sqrt(M) - 1)
    ps = zeros(N,L+1,K);
    for k = 1:K            
        for l = 0:L
            for n = 1:N
                #[l*l + 1, l*l + 2*l + 1]
                for i = 1:(2*l+1)                    
                    ps[n,l+1,k] += Ar[n,l*l + i,k]*Ar[n,l*l + i,k] + Ai[n,l*l + i,k]*Ai[n,l*l + i,k]                    
                end
            end
        end
    end

    if (dAr === nothing) | (dAi === nothing) 
        return ps 
    end

    dim, Nij, M, K = size(dAr)
    dps = zeros(dim,Nij,L+1,K);    
    for k = 1:K            
        for l = 0:L
            for n = 1:Nij
                ii = ai[n]
                ij = aj[n]
                for d = 1:dim
                    dtmp = 0.0
                    for i = 1:(2*l+1)                    
                        dtmp += dAr[d,n,l*l + i,k]*Ar[ii,l*l + i,k] + 
                                Ar[ii,l*l + i,k]*dAr[d,n,l*l + i,k] + 
                                dAi[d,n,l*l + i,k]*Ai[ii,l*l + i,k] +
                                Ai[ii,l*l + i,k]*dAi[d,n,l*l + i,k] 
                    end
                    dps[d,n,l+1,k] += dtmp
                end
            end
        end
    end

    return ps, dps   
end

function bispectrum1(Ar, Ai, dAr, dAi, cg, indl, indm, rowm, ai, aj)

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
    dbs = zeros(dim,Nij,J,K);    
    for k = 1:K            
        for j = 1:J
            l2 = indl[j,1];
            l1 = indl[j,2];
            l  = indl[j,3];             
            nm = rowm[j+1]-rowm[j];     
            for n = 1:Nij
                ii = ai[n]
                ij = aj[n]
                for d = 1:dim
                    dtmp = 0.0;
                    for i = 1:nm       
                        q = rowm[j]+i
                        c = cg[q];
                        m2 = indm[q,1];
                        m1 = indm[q,2];
                        m  = indm[q,3];
                        a1 = Ar[ii,m + l + (l*l) + 1,k]
                        b1 = Ai[ii,m + l + (l*l) + 1,k]
                        a2 = Ar[ii,m1 + l1 + (l1*l1) + 1,k]
                        b2 = Ai[ii,m1 + l1 + (l1*l1) + 1,k]
                        a3 = Ar[ii,m2 + l2 + (l2*l2) + 1,k]
                        b3 = Ai[ii,m2 + l2 + (l2*l2) + 1,k]

                        da1 = dAr[d,n,m + l + (l*l) + 1,k]
                        db1 = dAi[d,n,m + l + (l*l) + 1,k]
                        da2 = dAr[d,n,m1 + l1 + (l1*l1) + 1,k]
                        db2 = dAi[d,n,m1 + l1 + (l1*l1) + 1,k]
                        da3 = dAr[d,n,m2 + l2 + (l2*l2) + 1,k]
                        db3 = dAi[d,n,m2 + l2 + (l2*l2) + 1,k]

                        #tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);   
                        dtmp += c*(da1*a2*a3 + a1*da2*a3 + a1*a2*da3 +
                                        da2*b1*b3 + a2*db1*b3 + a2*b1*db3 + 
                                        da3*b1*b2 + a3*db1*b2 + a3*b1*db2 - 
                                        da1*b2*b3 - a1*db2*b3 - a1*b2*db3);       
                    end
                    dbs[d,n,j,k] += dtmp
                end
            end
        end
    end

    return bs, dbs   
end


function powerspectrum(A)

    N, M, K = size(A)
    L = Int32(sqrt(M) - 1)
    ps = zeros(N,L+1,K);
    for n = 1:N
        for k = 1:K            
            for l = 0:L
                #[l*l + 1, l*l + 2*l + 1]
                for i = 1:(2*l+1)                    
                    ps[n,l+1,k] += real(A[n,l*l + i,k])*real(A[n,l*l + i,k]) + imag(A[n,l*l + i,k])*imag(A[n,l*l + i,k])
                end
            end
        end
    end
    return ps 
end

function powerspectrum1(Ar, Ai)

    N, M, K = size(Ar)
    L = Int32(sqrt(M) - 1)
    ps = zeros(N,L+1,K);
    for n = 1:N
        for k = 1:K            
            for l = 0:L
                #[l*l + 1, l*l + 2*l + 1]
                for i = 1:(2*l+1)                    
                    ps[n,l+1,k] += Ar[n,l*l + i,k]*Ar[n,l*l + i,k] + Ai[n,l*l + i,k]*Ai[n,l*l + i,k]
                end
            end
        end
    end
    return ps 
end

function powerspectrum1(Ar, Ai, dAr, dAi)

    M, K = size(Ar)
    L = Int32(sqrt(M) - 1)
    ps = zeros(L+1,K);
    for k = 1:K            
        for l = 0:L
            #[l*l + 1, l*l + 2*l + 1]
            for i = 1:(2*l+1)                    
                ps[l+1,k] += Ar[l*l + i,k]*Ar[l*l + i,k] + Ai[l*l + i,k]*Ai[l*l + i,k]
            end
        end
    end

    if (dAr === nothing) | (dAi === nothing) 
        return ps 
    end

    dim, Nij, M, K = size(dAr)
    dps = zeros(dim,Nij,L+1,K);    
    for k = 1:K            
        for l = 0:L
            for n = 1:Nij
                for d = 1:dim
                    for i = 1:(2*l+1)                    
                        dps[d,n,l+1,k] += dAr[d,n,l*l + i,k]*Ar[l*l + i,k] + 
                                          Ar[l*l + i,k]*dAr[d,n,l*l + i,k] + 
                                          dAi[d,n,l*l + i,k]*Ai[l*l + i,k] +
                                          Ai[l*l + i,k]*dAi[d,n,l*l + i,k] 
                    end
                end
            end
        end
    end

    return ps, dps   
end

function powerspectrum2(Ar, Ai)

    N, M, K = size(Ar)
    L = Int32(sqrt(M) - 1)
    K2 = Int32(K*(K+1)/2)
    ps = zeros(N,L+1,K2);
    for n = 1:N
        kk = 0
        for k = 1:K     
            for k1 = k:K           
                kk = kk + 1 
                for l = 0:L
                    for i = 1:(2*l+1)                    
                        ps[n,l+1,kk] += Ar[n,l*l + i,k]*Ar[n,l*l + i,k1] + Ai[n,l*l + i,k]*Ai[n,l*l + i,k1]
                    end
                end
            end
        end
    end
    return ps 
end

function powerspectrum2(Ar, Ai, dAr, dAi)

    M, K = size(Ar)
    L = Int32(sqrt(M) - 1)
    K2 = Int32(K*(K+1)/2)
    ps = zeros(L+1,K2);
    kk = 0
    for k = 1:K     
        for k1 = k:K           
            kk = kk + 1 
            for l = 0:L
                for i = 1:(2*l+1)                    
                    ps[l+1,kk] += Ar[l*l + i,k]*Ar[l*l + i,k1] + Ai[l*l + i,k]*Ai[l*l + i,k1]
                end
            end
        end
    end

    if (dAr === nothing) | (dAi === nothing) 
        return ps 
    end

    dim, Nij, M, K = size(dAr)
    dps = zeros(dim,Nij,L+1,K2);    
    kk = 0
    for k = 1:K     
        for k1 = k:K         
            kk = kk + 1    
            for l = 0:L
                for n = 1:Nij
                    for d = 1:dim
                        for i = 1:(2*l+1)                    
                            dps[d,n,l+1,kk] += dAr[d,n,l*l + i,k]*Ar[l*l + i,k1] + 
                                            Ar[l*l + i,k]*dAr[d,n,l*l + i,k1] + 
                                            dAi[d,n,l*l + i,k]*Ai[l*l + i,k1] +
                                            Ai[l*l + i,k]*dAi[d,n,l*l + i,k1] 
                        end
                    end
                end
            end
        end
    end

    return ps, dps           
end

function bispectrum1(Ar, Ai, cg, indl, indm, rowm)

    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K);
    for n = 1:N
        for k = 1:K     
            for j = 1:J
                l2 = indl[j,1];
                l1 = indl[j,2];
                l  = indl[j,3];             
                tmp = 0;
                nm = rowm[j+1]-rowm[j];        
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
    return bs 
end

function bispectrum1(Ar, Ai, dAr, dAi, cg, indl, indm, rowm)

    J = size(indl,1);
    M, K = size(Ar)
    bs = zeros(J,K);
    for k = 1:K     
        for j = 1:J
            l2 = indl[j,1];
            l1 = indl[j,2];
            l  = indl[j,3];             
            tmp = 0;
            nm = rowm[j+1]-rowm[j];        
            for i = 1:nm
                q = rowm[j]+i
                c = cg[q];
                m2 = indm[q,1];
                m1 = indm[q,2];
                m  = indm[q,3];
                a1 = Ar[m + l + (l*l) + 1,k]
                b1 = Ai[m + l + (l*l) + 1,k]
                a2 = Ar[m1 + l1 + (l1*l1) + 1,k]
                b2 = Ai[m1 + l1 + (l1*l1) + 1,k]
                a3 = Ar[m2 + l2 + (l2*l2) + 1,k]
                b3 = Ai[m2 + l2 + (l2*l2) + 1,k]
                tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                         
            end                                 
            bs[j,k] = tmp
        end
    end

    if (dAr === nothing) | (dAi === nothing) 
        return bs 
    end

    dim, Nij, M, K = size(dAr)
    dbs = zeros(dim,Nij,J,K);    
    for k = 1:K            
        for j = 1:J
            l2 = indl[j,1];
            l1 = indl[j,2];
            l  = indl[j,3];             
            nm = rowm[j+1]-rowm[j];     
            for n = 1:Nij
                for d = 1:dim
                    for i = 1:nm       
                        q = rowm[j]+i
                        c = cg[q];
                        m2 = indm[q,1];
                        m1 = indm[q,2];
                        m  = indm[q,3];
                        a1 = Ar[m + l + (l*l) + 1,k]
                        b1 = Ai[m + l + (l*l) + 1,k]
                        a2 = Ar[m1 + l1 + (l1*l1) + 1,k]
                        b2 = Ai[m1 + l1 + (l1*l1) + 1,k]
                        a3 = Ar[m2 + l2 + (l2*l2) + 1,k]
                        b3 = Ai[m2 + l2 + (l2*l2) + 1,k]

                        da1 = dAr[d,n,m + l + (l*l) + 1,k]
                        db1 = dAi[d,n,m + l + (l*l) + 1,k]
                        da2 = dAr[d,n,m1 + l1 + (l1*l1) + 1,k]
                        db2 = dAi[d,n,m1 + l1 + (l1*l1) + 1,k]
                        da3 = dAr[d,n,m2 + l2 + (l2*l2) + 1,k]
                        db3 = dAi[d,n,m2 + l2 + (l2*l2) + 1,k]

                        #tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);   
                        dbs[d,n,j,k] += c*(da1*a2*a3 + a1*da2*a3 + a1*a2*da3 +
                                        da2*b1*b3 + a2*db1*b3 + a2*b1*db3 + 
                                        da3*b1*b2 + a3*db1*b2 + a3*b1*db2 - 
                                        da1*b2*b3 - a1*db2*b3 - a1*b2*db3);       
                    end
                end
            end
        end
    end

    return bs, dbs   
end

function bispectrum2(Ar, Ai, cg, indl, indm, rowm)

    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K*K);
    for n = 1:N
        for k = 1:K   
            for k1 = 1:K     
                for j = 1:J
                    l2 = indl[j,1];
                    l1 = indl[j,2];
                    l  = indl[j,3];             
                    tmp = 0;
                    nm = rowm[j+1]-rowm[j];        
                    for i = 1:nm
                        q = rowm[j]+i
                        c = cg[q];
                        m2 = indm[q,1];
                        m1 = indm[q,2];
                        m  = indm[q,3];
                        a1 = Ar[n,m + l + (l*l) + 1,k]
                        b1 = Ai[n,m + l + (l*l) + 1,k]
                        a2 = Ar[n,m1 + l1 + (l1*l1) + 1,k1]
                        b2 = Ai[n,m1 + l1 + (l1*l1) + 1,k1]
                        a3 = Ar[n,m2 + l2 + (l2*l2) + 1,k1]
                        b3 = Ai[n,m2 + l2 + (l2*l2) + 1,k1]
                        tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                         
                    end                                 
                    bs[n,j,k1 + (k-1)*K] = tmp
                end
            end
        end
    end
    return bs 
end


function bispectrum3(Ar, Ai, cg, indl, indm, rowm)

    J = size(indl,1);
    N, M, K = size(Ar)
    bs = zeros(N,J,K*K*K);
    for n = 1:N
        for k = 1:K   
            for k1 = 1:K 
                for k2 = 1:K     
                    for j = 1:J
                        l2 = indl[j,1];
                        l1 = indl[j,2];
                        l  = indl[j,3];             
                        tmp = 0;
                        nm = rowm[j+1]-rowm[j];        
                        for i = 1:nm
                            q = rowm[j]+i
                            c = cg[q];
                            m2 = indm[q,1];
                            m1 = indm[q,2];
                            m  = indm[q,3];
                            a1 = Ar[n,m + l + (l*l) + 1,k]
                            b1 = Ai[n,m + l + (l*l) + 1,k]
                            a2 = Ar[n,m1 + l1 + (l1*l1) + 1,k1]
                            b2 = Ai[n,m1 + l1 + (l1*l1) + 1,k1]
                            a3 = Ar[n,m2 + l2 + (l2*l2) + 1,k2]
                            b3 = Ai[n,m2 + l2 + (l2*l2) + 1,k2]
                            tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                         
                        end                                 
                        bs[n,j,k2 + (k1-1)*K + (k-1)*K*K] = tmp
                    end
                end
            end
        end
    end
    return bs 
end




