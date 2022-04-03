function makejk(n)
    
    m = Int32((n-1)*n/2);
    indj = zeros(Int32,m);
    indk = zeros(Int32,m);
    k1 = 1;
    for i = 1:(n-1)
        ind = k1:(k1+(n-i)-1);
        indj[ind] .= i;
        indk[ind] = (i+1):n;
        k1 = k1 + (n-i);
    end
    
    return indj, indk
        
    end
    
function makejk!(indj, indk, n)
    
    k1 = 1;
    for i = 1:(n-1)
        ind = k1:(k1+(n-i)-1);
        indj[ind] .= i;
        indk[ind] = (i+1):n;
        k1 = k1 + (n-i);
    end
    
    #return indj, indk        
end
    

function makejk(e2ij, f2ij, pairnum, tripletnum, ilist)

    inum = length(ilist);
    ijknum = sum(tripletnum[ilist.+1]-tripletnum[ilist]);        
    M = size(e2ij,2)
    uij = zeros(ijknum,M);
    uik = zeros(ijknum,M);      
    wij = zeros(3,ijknum,M);
    wik = zeros(3,ijknum,M);      
    n = maximum(pairnum[ilist.+1]-pairnum[ilist]);  
    m = Int32((n-1)*n/2);      
    indj = zeros(Int32, m)
    indk = zeros(Int32, m)
    ind = zeros(Int32, m)
    ind2 = zeros(Int32, n)    
    ei = zeros(n, M)
    f1 = zeros(n, M)
    f2 = zeros(n, M)
    f3 = zeros(n, M)
    count = 0;
    count2 = 0
    for ii=1:inum
        i = ilist[ii]; # atom i        
        n2 = pairnum[i+1] - pairnum[i]; 
        #ind2[1:n2] = (count2+1):(count2+n2);                
        for i2 = 1:n2
            ind2[i2] = count2+i2 
        end
        #indj, indk = makejk(n2);        
        makejk!(indj, indk, n2);
        #m = Int32((n2-1)*n2/2);      

        n = tripletnum[i+1] - tripletnum[i];         
        #ind = (count+1):(count+n);      
        for i2 = 1:n
            ind[i2] = count+i2 
        end

        for i1 = 1:M
            for i2 = 1:n2
                ei[i2,i1] = e2ij[ind2[i2],i1]
                f1[i2,i1] = f2ij[1+3*(ind2[i2]-1),i1]
                f2[i2,i1] = f2ij[2+3*(ind2[i2]-1),i1]
                f3[i2,i1] = f2ij[3+3*(ind2[i2]-1),i1]
            end
        end
        for i1 = 1:M
            for i2 = 1:n
                uij[ind[i2],i1] = ei[indj[i2],i1]
                uik[ind[i2],i1] = ei[indk[i2],i1]
                wij[1,ind[i2],i1] = f1[indj[i2],i1]
                wik[1,ind[i2],i1] = f1[indk[i2],i1]       
                wij[2,ind[i2],i1] = f2[indj[i2],i1]
                wik[2,ind[i2],i1] = f2[indk[i2],i1]       
                wij[3,ind[i2],i1] = f3[indj[i2],i1]
                wik[3,ind[i2],i1] = f3[indk[i2],i1]         
            end
        end

        # uij[ind,:] = ei[indj[1:m],:]
        # uik[ind,:] = ei[indk[1:m],:]
        # wij[1,ind,:] = f1[indj[1:m],:]
        # wik[1,ind,:] = f1[indk[1:m],:]
        # wij[2,ind,:] = f2[indj[1:m],:]
        # wik[2,ind,:] = f2[indk[1:m],:]
        # wij[3,ind,:] = f3[indj[1:m],:]
        # wik[3,ind,:] = f3[indk[1:m],:]
        
        # #ei = e2ij[ind2,:]
        # ei[1:n2,:] = e2ij[ind2,:]
        # uij[ind,:] = ei[indj[1:m],:]
        # uik[ind,:] = ei[indk[1:m],:]
        
        # #fi = f2ij[1,ind2,:]
        # ei[1:n2,:] = f2ij[1,ind2,:]
        # wij[1,ind,:] = ei[indj[1:m],:]
        # wik[1,ind,:] = ei[indk[1:m],:]
        # #fi = f2ij[2,ind2,:]
        # ei[1:n2,:] = f2ij[2,ind2,:]
        # wij[2,ind,:] = ei[indj[1:m],:]
        # wik[2,ind,:] = ei[indk[1:m],:]
        # #fi = f2ij[3,ind2,:]
        # ei[1:n2,:] = f2ij[3,ind2,:]
        # wij[3,ind,:] = ei[indj[1:m],:]
        # wik[3,ind,:] = ei[indk[1:m],:]

        count = count + n;                   
        count2 = count2 + n2;                        
    end
    
    return uij, uik, wij, wik      
end
        
function neightripletlist(pairlist, pairnumsum, ilist)

    pairnum = Int32.(pairnumsum[2:end]-pairnumsum[1:end-1]);
    tripletnum  = Int32.((pairnum.-1).*pairnum/2);

    inum = length(ilist)
    nt = 0
    for i = 1:inum
        n = pairnum[i]
        nt = nt + (n-1)*n/2
    end
    nt = Int32(nt)
    tripletlistj = zeros(Int32,nt);
    tripletlistk = zeros(Int32,nt);

    nt = 0
    for ii = 1:inum
        n1 = pairnumsum[ii];
        n2 = pairnumsum[ii+1];        
        g = pairlist[(n1+1):n2]; # a list of neighbors around atom i
        
        n = pairnum[ii]
        indj, indk = makejk(n);        
        # tripletlistj = [tripletlistj; g[indj]];      
        # tripletlistk = [tripletlistk; g[indk]];     
        k = Int32((n-1)*n/2)  
        tripletlistj[(nt+1):(nt+k)] =  g[indj];      
        tripletlistk[(nt+1):(nt+k)] =  g[indk];      
        nt = nt + k
    end

    tripletlist = [tripletlistj[:] tripletlistk[:]];
    tripletnum = [0; cumsum(tripletnum)];

    tripletlist = Int32.(tripletlist)
    tripletnum  = Int32.(tripletnum)

    return tripletlist, tripletnum

end
    
    