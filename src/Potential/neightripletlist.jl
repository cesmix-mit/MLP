function makejk(n)
    
m = Int64((n-1)*n/2);
indj = zeros(Int64,m);
indk = zeros(Int64,m);
k1 = 1;
for i = 1:(n-1)
    ind = k1:(k1+(n-i)-1);
    indj[ind] .= i;
    indk[ind] = (i+1):n;
    k1 = k1 + (n-i);
end

return indj, indk
    
end

function makejk(e2ij, f2ij, pairnum, tripletnum, ilist)

    inum = length(ilist);
    ijknum = sum(tripletnum[ilist.+1]-tripletnum[ilist]);        
    M = size(e2ij,2)
    uij = zeros(ijknum,M);
    uik = zeros(ijknum,M);      
    wij = zeros(3,ijknum,M);
    wik = zeros(3,ijknum,M);      
    count = 0;
    count2 = 0
    for ii=1:inum
        i = ilist[ii]; # atom i        
        n2 = pairnum[i+1] - pairnum[i]; 
        ind2 = (count2+1):(count2+n2);                
        indj, indk = makejk(n2);        

        n = tripletnum[i+1] - tripletnum[i];         
        ind = (count+1):(count+n);      
        
        ei = e2ij[ind2,:]
        uij[ind,:] = ei[indj,:]
        uik[ind,:] = ei[indk,:]
        
        fi = f2ij[1,ind2,:]
        wij[1,ind,:] = fi[indj,:]
        wik[1,ind,:] = fi[indk,:]
        fi = f2ij[2,ind2,:]
        wij[2,ind,:] = fi[indj,:]
        wik[2,ind,:] = fi[indk,:]
        fi = f2ij[3,ind2,:]
        wij[3,ind,:] = fi[indj,:]
        wik[3,ind,:] = fi[indk,:]

        count = count + n;                   
        count2 = count2 + n2;                        
    end
    
    return uij, uik, wij, wik      
end
        
function neightripletlist(pairlist, pairnumsum, ilist)

pairnum = Int64.(pairnumsum[2:end]-pairnumsum[1:end-1]);
tripletnum  = Int64.((pairnum.-1).*pairnum/2);

inum = length(ilist)
nt = 0
for i = 1:inum
    n = pairnum[i]
    nt = nt + (n-1)*n/2
end
nt = Int64(nt)
tripletlistj = zeros(Int64,nt);
tripletlistk = zeros(Int64,nt);

nt = 0
for ii = 1:inum
    n1 = pairnumsum[ii];
    n2 = pairnumsum[ii+1];        
    g = pairlist[(n1+1):n2]; # a list of neighbors around atom i
    
    n = pairnum[ii]
    indj, indk = makejk(n);        
    # tripletlistj = [tripletlistj; g[indj]];      
    # tripletlistk = [tripletlistk; g[indk]];     
    k = Int64((n-1)*n/2)  
    tripletlistj[(nt+1):(nt+k)] =  g[indj];      
    tripletlistk[(nt+1):(nt+k)] =  g[indk];      
    nt = nt + k
end

tripletlist = [tripletlistj[:] tripletlistk[:]];
tripletnum = [0; cumsum(tripletnum)];

tripletlist = Int64.(tripletlist)
tripletnum  = Int64.(tripletnum)

return tripletlist, tripletnum

end
    
    