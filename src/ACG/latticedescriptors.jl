using Preprocessing, Potential, LinearAlgebra

function PODdescriptors(species = [:X], rin = 1.0, rcut = 5.0, gamma = [0.0, 2, 4])
    
descriptors[1] = Potential.POD(nbody=1, species = species, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = Potential.POD(nbody=2, species = species, pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = Potential.POD(nbody=3, species = species, pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors)

return descriptors

end

function latticedescriptors(a1, a2, a3, atombasis, atomtype, descriptors)

    pbc = [1, 1, 1]
    K = size(a1,2)    
    N = size(atombasis,3)    
    d, fatom, eatom = Potential.podefatom(atombasis[:,:,1], atomtype, a1[:,1], a2[:,1], a3[:,1], pbc, descriptors)
    natom, M = size(eatom)    
    gdesc = zeros(M, N)   
    edesc = zeros(natom*M, N)    
    fdesc = zeros(3*natom*M, N)    
    gdesc[:,1] = d[:]    
    edesc[:,1] = eatom[:] 
    fdesc[:,1] = fatom[:]     
    for n = 2:N
        display(n)
        k = (K == 1) ? 1 : n;
        d, fatom, eatom = Potential.podefatom(atombasis[:,:,n], atomtype, a1[:,k], a2[:,k], a3[:,k], pbc, descriptors)
        gdesc[:,n] = d[:]    
        edesc[:,n] = eatom[:] 
        fdesc[:,n] = fatom[:] 
    end

    return edesc, fdesc, gdesc   
end

function latticepoddescriptors(a1, a2, a3, atombasis, atomtype, descriptors)

    pbc = [1, 1, 1]
    K = size(a1,2)    
    N = size(atombasis,3)    
    d, eatom = Potential.podeatom(atombasis[:,:,1], atomtype, a1[:,1], a2[:,1], a3[:,1], pbc, descriptors)
    natom, M = size(eatom)    
    gdesc = zeros(M, N)   
    edesc = zeros(natom*M, N)    
    gdesc[:,1] = d[:]    
    edesc[:,1] = eatom[:] 
    for n = 2:N
        display(n)
        k = (K == 1) ? 1 : n;
        d, eatom = Potential.podeatom(atombasis[:,:,n], atomtype, a1[:,k], a2[:,k], a3[:,k], pbc, descriptors)
        gdesc[:,n] = d[:]    
        edesc[:,n] = eatom[:] 
    end

    return edesc, gdesc   
end


function latticesimilaritymatrices(a1, a2, a3, atombasis, atomtype, descriptors)

    edesc, fdesc, gdesc  = latticedescriptors(a1, a2, a3, atombasis, atomtype, descriptors)

    E = similaritymeasure(edesc)    
    F = similaritymeasure(fdesc)    
    G = similaritymeasure(gdesc)    
    # F = 0.0
    # G = 0.0

    return E, F, G, edesc, fdesc, gdesc  
end


function similaritymatrices(a1, a2, a3, atombasis, atomtype, descriptors)

    pbc = [1, 1, 1]
    K = size(a1,2)    
    N = size(atombasis,3)    
    d, eatom = Potential.podeatom(atombasis[:,:,1], atomtype, a1[:,1], a2[:,1], a3[:,1], pbc, descriptors)
    natom, M = size(eatom)    
    gdesc = zeros(M, N)   
    edesc = zeros(natom*M, N)    
    gdesc[:,1] = d[:]    
    edesc[:,1] = eatom[:] 
    for n = 2:N
        display(n)
        k = (K == 1) ? 1 : n;
        d, eatom = Potential.podeatom(atombasis[:,:,n], atomtype, a1[:,k], a2[:,k], a3[:,k], pbc, descriptors)
        gdesc[:,n] = d[:]    
        edesc[:,n] = eatom[:] 
    end

    E = similaritymeasure(edesc)    
    G = similaritymeasure(gdesc)    

    return E, G, edesc, gdesc   
end

function distancemeasure(M)
    return [ norm(a-b) for a in eachcol(M), b in eachcol(M) ]
end

function gaussiandistancemeasure(M)
    return [ exp(-norm(a/norm(a)-b/norm(b))^2) for a in eachcol(M), b in eachcol(M) ]
end

function cosmeasure(M)    
    return [ sum(a.*b)/(norm(a)*norm(b)) for a in eachcol(M), b in eachcol(M) ]  
end

function sinmeasure(M)
    CosA = [ sum(a.*b)/(norm(a)*norm(b)) for a in eachcol(M), b in eachcol(M) ]    
    return sqrt.( 1.0 .- CosA.^2)
end

function similaritymeasure(M)    
    l, n = size(M)
    for i = 1:n
        M[:,i] = M[:,i]/norm(M[:,i])
    end    
    
    # G = [ exp(-norm(a-b)^2) for a in eachcol(M), b in eachcol(M) ]
    # CosA = [abs(sum(a.*b) for a in eachcol(M), b in eachcol(M) ]        
    # S = G.*CosA 
    # S = zeros(n,n)
    # for i = 1:n
    #     display(i)
    #     for j = i:n         
    #         ta = 0.0;
    #         #tb = 0.0;
    #         for k = 1:l                
    #             ta += M[k,i]*M[k,j]; 
    #             #tb += (M[k,i]-M[k,j])*(M[k,i]-M[k,j]);
    #         end
    #         tm = abs(ta);#*exp(-tb);            
    #         S[i,j] = tm;
    #         S[j,i] = tm; 
    #     end
    # end

    # S2 = M'*M;
    # display(S)
    # display(S2)
    # display(maximum(abs.(S[:]-S2[:])))
    # error("here")

    S = M'*M;
    return S
end


