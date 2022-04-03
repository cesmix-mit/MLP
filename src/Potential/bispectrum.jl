function radialsphericalharmonicbasis_sum(phir, phii, ai, N)

    Nij, M, K = size(phir)
    Ar = zeros(N, M, K)
    Ai = zeros(N, M, K)
    for k = 1:K     
        for m = 1:M
            for n = 1:Nij 
                ii = ai[n]
                Ar[ii,m,k] += phir[n,m,k]
                Ai[ii,m,k] += phii[n,m,k]
            end
        end
    end              

    return Ar, Ai
end

function powerspectrum(Ar, Ai, dAr, dAi, ai, aj)

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
    dps = zeros(dim,N,L+1,K);    
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
                    dps[d,ii,l+1,k] += dtmp
                    dps[d,ij,l+1,k] -= dtmp 
                end
            end
        end
    end

    return ps, dps   
end

function bispectrum(Ar, Ai, dAr, dAi, cg, indl, indm, rowm, ai, aj)

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
    dbs = zeros(dim,N,J,K);        
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
                    dbs[d,ii,j,k] += dtmp
                    dbs[d,ij,j,k] -= dtmp 
                end
            end
        end
    end

    return bs, dbs   
end
