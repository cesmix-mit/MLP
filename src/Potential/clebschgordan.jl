function cgconditions(L)

    inc = 0;
    for l1 = 0:L
        for l2 = 0:L
            for l = 0:L
                if  (abs(l1-l2) <= l) & (l <= l1 + l2)  & (abs(l-l2) <= l1) & 
                    (l1 <= l + l2)  & (abs(l-l1) <= l2) & (l2 <= l + l1) & 
                    (l2>=l1) & (l1>=l) & (rem(l1+l2+l,2)==0)           
                    inc = inc + 1;
                end
            end
        end
    end
    
    ind = zeros(Int64,inc,3);
    inc = 0;
    for l1 = 0:L
        for l2 = 0:L
            for l = 0:L
                if  (abs(l1-l2) <= l) & (l <= l1 + l2) & (abs(l-l2) <= l1) & 
                    (l1 <= l + l2) & (abs(l-l1) <= l2) & (l2 <= l + l1) & 
                    (l2>=l1) & (l1>=l) & (rem(l1+l2+l,2)==0)           
                        inc = inc + 1;
                        ind[inc,:] = [l2, l1, l];             
                end
            end
        end
    end

    return ind
end

function fac(n)    
    return factorial(big(n))
end

function fac(n::Vector{Float64})    
    m = length(n)
    f = zeros(m)
    for i = 1:m
        f[i] = factorial(big(n[i]))
    end
    return f
end

function clebschgordan(j1,m1,j2,m2,j,m)
    
    if (m1+m2-m != 0) | (j < abs(j1-j2)) | (j > j1+j2)
        C = 0;
        return    
    end
        
    #---compute valid k values for summation---%
    k1 = maximum([0,j2-j-m1,j1-j+m2])
    k2 = minimum([j1+j2-j,j1-m1,j2+m2])
    k = Float64.(Array(k1:k2))
    
    #---compute coefficient---%
    C = sqrt((2*j+1)*fac(j+j1-j2)*fac(j+j2-j1)*fac(j1+j2-j)/fac(j1+j2+j+1))*
        sqrt(fac(j+m)*fac(j-m)*fac(j1+m1)*fac(j1-m1)*fac(j2+m2)*fac(j2-m2))*
        sum(((-1).^k)./(fac(k).*fac(j1+j2-j.-k).*fac(j1-m1.-k).*fac(j2+m2.-k).*fac(j-j2+m1.+k).*fac(j-j1-m2.+k)));
        
    return C;
end

function cgcoefficients(ind)
    
    N = size(ind,1);
    
    nc = 0;
    for n = 1:N
        l2 = ind[n,1];
        l1 = ind[n,2];
        l = ind[n,3];
        for m = -l:1:l
            for m1 = -l1:1:l1
                for m2 = -l2:1:l2
                    if (m1+m2-m == 0)  
                        cg = clebschgordan(l1,m1,l2,m2,l,m);
                        if abs(cg)>1e-12
                            nc = nc + 1;
                        end
                    end
                end
            end
        end
    end
    
    c = zeros(nc);
    indm = zeros(Int64,nc,3);
    nc = 1;
    rowm = zeros(Int64,N);
    for n = 1:N
        l2 = ind[n,1];
        l1 = ind[n,2];
        l = ind[n,3];
        count = 0;
        for m = -l:1:l
            for m1 = -l1:1:l1
                for m2 = -l2:1:l2                    
                    if (m1+m2-m == 0)  
                        cg = clebschgordan(l1,m1,l2,m2,l,m);
                        if abs(cg)>1e-12
                            c[nc] = cg;
                            indm[nc,:] = [m2, m1, m];   
                            nc = nc + 1;
                            count = count + 1;                        
                        end
                    end                        
                end
            end
        end
        rowm[n] = count;
    end
    rowm = cumsum([0; rowm]);

    return c,indm,rowm
end

# L = 4;
# indl = cgconditions(L);
# cg,indm,rowm = cgcoefficients(indl);

