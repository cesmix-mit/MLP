function numlattices(rcutmax, pbc, a, b, c)

    d = (length(c) < 3) ? 2 : 3

    m = 0
    n = 0
    p = 0

    if pbc[1] == 1
        m = Int32(ceil(rcutmax/a[1]))
    end
    if pbc[2] == 1
        n = Int32(ceil(rcutmax/b[2]))
    end
    if d==3
        if pbc[3] == 1
            p = Int32(ceil(rcutmax/c[3]))
        end
    end
    # m = m + 1;
    # n = n + 1;
    # p = p + 1;

    return m, n, p
end

function referenceorigins(m, n, p, d)

xref = zeros(d, (2*m+1)*(2*n+1)*(2*p+1))
for i = 1:(2*p+1) # z
    for j = 1:(2*n+1) # y
        for k = 1:(2*m+1) # x
            t = k + (2*m+1)*(j-1) + (2*m+1)*(2*n+1)*(i-1)
            xref[1, t] = k-m-1;  
            xref[2, t] = j-n-1;  
            if d==3 
                xref[3, t] = i-p-1;  
            end
        end
    end
end

return xref

end

function latticeorigins(a, b, c, m, n, p)

    d = (length(c) < 3) ? 2 : 3
 
    xo = referenceorigins(m, n, p, d)

    if d==2
        B2R, R2B = domainmaps(a, b);    
    else
        B2R, R2B = domainmaps(a, b, c);    
    end

    latorigins = R2B*xo; 

    return latorigins
end

function latticecoords(x, rcutmax, pbc, a, b, c)

    m, n, p = numlattices(rcutmax, pbc, a, b, c)
    
    latorigins = latticeorigins(a, b, c, m, n, p)
    ind = (m+1) + (2*m+1)*(n) + (2*m+1)*(2*n+1)*(p)

    # display([m n p])
    # display(latorigins)
    # error("here")

    d = size(x,1)
    nx = size(x,2)
    nl = size(latorigins,2)  # number of lattices

    y = zeros(d, nx*nl)
    for j = 1:nx
        for k = 1:d
            y[k,j] = x[k,j] 
        end
    end
    q = nx;
    for i = 1:nl
        if i != ind            
            for j = 1:nx
                q = q + 1
                for k = 1:d
                    y[k,q] = latorigins[k,i] + x[k,j]
                end
            end
        end
    end

    alist = zeros(Int32,nx*nl)
    for i = 1:nl
        for j = 1:nx
            alist[j + nx*(i-1)] = j
        end
    end

    # # get 
    # neighi, neighj, neighnum = neighborlist(y, rcutmax);

    # neighi = neighi[neighi .<= nx];
    # ne = length(neighi);

    # # list of neighbors for each atom i in x
    # neighlist = neighj[1:ne];

    # # number of neighbors for each atom i in x
    # neighnum = neighnum[1:nx];
    # neighnum = [0; cumsum(neighnum)];

    neighi, neighlist, neighnum = neighborlist(y, rcutmax, nx)
    neighnumsum = Int32.([0; cumsum(neighnum)]);
    
    #display(neighnum)
    # display(neighlist)
    # display(y)
    # display(a)
    # display(b)
    # display(c)
    # display(pbc)

# int fullneighborlist(double *y, int *alist, int *selflist, int *neighlist, int *numneigh, 
#         int *numneighsum, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);
    # y1 = zeros(d, nx*nl)
    # alist1 = zeros(Int32,nx*nl);
    # neighlist1 = zeros(Int32,nx*nx*nl);
    # numneigh = zeros(Int32,nx);
    # numneighsum = zeros(Int32,nx+1);
    # count = podfullneighborlist(y1, alist1, neighlist1, numneigh, numneighsum, x, a, b, c, rcutmax, Int32.(pbc), Int32(nx));

    # k = length(neighlist)    
    # display([alist[:]-alist1[:].-1])
    # display(y-y1)    
    # display(neighlist-neighlist1[1:k].-1)    
    # display(neighnum-numneigh)
    # display([neighnumsum[:]-numneighsum[:]])
    # display([k count])

    return y, alist, neighlist, neighnumsum, neighnum
end

