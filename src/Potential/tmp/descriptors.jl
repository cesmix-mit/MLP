using Preprocessing
using ACEpot

function ACEdescriptors(x, t, a, b, c, pbc, mldescriptors::Int64)
    return nothing, nothing, nothing
end

function mlglobaldescriptors(x, t, a, b, c, pbc, mldescriptors)
    if mldescriptors.name == "SNAP" 
        if typeof(mldescriptors) == Potential.SnaStruct
            rcutmax = 0.0
            for ielem = 1:mldescriptors.ntypes
                rcutmax = max(2.0*mldescriptors.radelem[1+ielem]*mldescriptors.rcutfac,rcutmax)
            end        
            mldescriptors.rcutmax = rcutmax 
        end
        mld, mldd, mldv = snapdescriptors(x, t, a, b, c, pbc, mldescriptors)
        sz = size(mldd)
        mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
    end
    if mldescriptors.name == "ACE" 
        mld, mldd, mldv = ACEpot.ACEdescriptors(x, t, a, b, c, pbc, mldescriptors)
        sz = size(mldd)
        mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
    end
    if mldescriptors.name == "POD"         
        mld, mldd, mldv = PODdescriptors(x, t, a, b, c, pbc, mldescriptors.rcut, mldescriptors)
        sz = size(mldd)        
        mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
    end
    return mld, mldd, mldv
end

function machinelearningdescriptors(x, t, a, b, c, pbc, mldescriptors)
    
    # bispectrum, bispectrum derivatives, bispectrum virial
    if length(mldescriptors) > 0
        mld = 0; mldd = 0; mldv = 0;
        for i = 1:length(mldescriptors)       
            md, mdd, mdv = mlglobaldescriptors(x, t, a, b, c, pbc, mldescriptors[i])   
            if i == 1
                mld = 1.0*md
                mldd = 1.0*mdd
                mldv = 1.0*mdv
            else                
                mld = [mld[:]; md[:]]
                mldd = cat(mldd, mdd, dims=2);    
                mldv = cat(mldv, mdv, dims=2);    
            end
        end
        return mld, mldd, mldv 
    else
        return nothing, nothing, nothing 
    end
end
    
function globaldescriptors(x, q, t, a, b, c, pbc, eta, kappa, emdescriptors, mldescriptors)

    if length(emdescriptors) > 0
        rcutmax = 0.0
        for i = 1:length(emdescriptors)                
            rcutmax = max(rcutmax, emdescriptors[i].rcut)        
        end
        ed, edd, edv = empiricaldescriptors(x, q, t, a, b, c, pbc, rcutmax, eta, kappa, emdescriptors)
        sz = size(edd)
        if length(sz)==3
            edd = reshape(edd, (sz[1]*sz[2], sz[3]))
        else
            edd = reshape(edd, (sz[1]*sz[2], 1))
        end
    end

    # bispectrum, bispectrum derivatives, bispectrum virial
    if length(mldescriptors) > 0
        mld, mldd, mldv = machinelearningdescriptors(x, t, a, b, c, pbc, mldescriptors)      
    end

    if (length(emdescriptors) > 0) & (length(mldescriptors) > 0)
        # ed = 10.0*ed 
        # edd = 4000.0*edd; 
        d = [ed[:]; mld[:]]
        dd = cat(edd, mldd, dims=2);    
        # ddnorm = sqrt.(sum(dd.*dd,dims=1))
        # show(stdout, "text/plain", d)
        # show(stdout, "text/plain", ddnorm)
        #display(d') 
        # display(ddnorm) 
        dv = cat(edv, mldv, dims=2);    
        return d, dd, dv    
    elseif length(emdescriptors) > 0
        return ed, edd, edv
    elseif length(mldescriptors) > 0
        return mld, mldd, mldv 
    else
        error("empirical descriptors and ML descriptors are empty")
    end    
end

function globaldescriptors(config, indices, eta, kappa, normalizeenergy, 
    normalizestress, descriptors, potential=nothing)

# getdescriptors
emdescriptors, mldescriptors = getemldescriptors(descriptors)


# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

nf = zeros(Int64, length(indices))
be = []
bf = []
bs  = []
Ae = []
Af = []
As  = []
for i = 1:length(indices)
    ci = indices[i]
    if ci > nconfigs
        error("Indices are greater than nconfigs")
    end

    x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
    t = config.t[:,(natom[ci]+1):natom[ci+1]] # atom types of configuration ci
    if length(config.q)>0
        q = config.q[:,(natom[ci]+1):natom[ci+1]]
    else
        q = []
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

    d, dd, dv = globaldescriptors(x, q, t, a[:], b[:], c[:], pbc[:], eta, kappa, emdescriptors, mldescriptors)        
    
    if normalizestress==1
        vol = latticevolume(a[:], b[:], c[:])
        dv = (1.6021765e6/vol) * dv
    end

    if (potential !== nothing)  
        potential = Preprocessing.deletenothing(potential)
        if (length(potential) > 0)             
            rcutmax = 0.0
            for j = 1:length(potential)                
                rcutmax = max(rcutmax, potential[j].rcut)            
            end
            eref, e, fref, v, sref = empiricalpotential(x, q, t, a[:], b[:], c[:], pbc[:], rcutmax, eta, kappa, potential)        
            if normalizestress==1
                vol = latticevolume(a[:], b[:], c[:])
                sref = (1.6021765e6/vol) * sref
            end
        else
            eref = 0.0; fref = [0.0]; sref = [0.0];
        end
    else
        eref = 0.0; fref = [0.0]; sref = [0.0];
    end

    if normalizeenergy==1
        d = reshape(d, (1,length(d)))/config.natom[ci]
    else
        d = reshape(d, (1,length(d)))
    end

    nf[i] = size(dd,1)
    if i == 1
        Ae = 1.0*d
        Af = 1.0*dd
        As = 1.0*dv
    else
        Ae = cat(Ae, d, dims=1)
        Af = cat(Af, dd, dims=1)
        As = cat(As, dv, dims=1)
    end

    if length(config.e)>0
        e = config.e[indices[i]] - eref
        if normalizeenergy==1
            e = e/config.natom[indices[i]]
        end
        be = [be; e]
    end    
    if length(config.f)>0
        f = config.f[:,(natom[ci]+1):natom[ci+1]]
        bf = [bf; f[:] .- fref[:]]    
    end
    if length(config.stress)>0
        if size(config.stress,1)==9
            ind = [1, 5, 9, 6, 3, 2]        
            s = config.stress[ind,ci]   
        else
            s = config.stress[:,ci]   
        end             
        bs = [bs; s[:] .- sref[:]]    
    end    
end

return Ae, Af, As, be, bf, bs, nf

end

function globaldescriptors(data, descriptors, options, potential=nothing)
        
    data = Preprocessing.deletenothing(data)
    descriptors = Preprocessing.deletenothing(descriptors)

    normalizeenergy = options.normalizeenergy
    normalizestress = options.normalizestress
    eta = options.eta 
    kappa = options.kappa 
    pbc = options.pbc 
    a = options.a
    b = options.b
    c = options.c
    
    if pbc === nothing        
        pbc = [1; 1; 1] # periodic boundary conditions
    end
    
    n = length(data)
    Ae = 1.0; Af = 1.0; As = 1.0; be = 1.0; bf = 1.0; bs = 1.0;
    nf = Array{Any}(nothing, n);
    for i = 1:n
        # training data set
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        training = Array(1:config.nconfigs)      

        Aei, Afi, Asi, bei, bfi, bsi, nfi = globaldescriptors(config, training, eta, 
                kappa, normalizeenergy, normalizestress, descriptors, potential)
        
        nf[i] = nfi 
        if (i == 1) 
            Ae = 1.0*Aei
            Af = 1.0*Afi
            As = 1.0*Asi
            be = 1.0*bei
            bf = 1.0*bfi
            bs = 1.0*bsi
        else         
            Ae, be = catarray(Ae, be, Aei, bei)
            Af, bf = catarray(Af, bf, Afi, bfi)
            As, bs = catarray(As, bs, Asi, bsi)                        
        end                        
    end

    return Ae, Af, As, be, bf, bs, nf
end

function globaldescriptors(data::Preprocessing.DataStruct, descriptors, options, potential=nothing)
        
    descriptors = Preprocessing.deletenothing(descriptors)

    normalizeenergy = options.normalizeenergy
    normalizestress = options.normalizestress
    eta = options.eta 
    kappa = options.kappa 
    pbc = options.pbc 
    a = options.a
    b = options.b
    c = options.c
    
    if pbc === nothing        
        pbc = [1; 1; 1] # periodic boundary conditions
    end
    
    config, indices = Preprocessing.readconfigdata(data, pbc, a, b, c)                  
    training = Array(1:config.nconfigs)      

    Ae, Af, As, be, bf, bs, nf = globaldescriptors(config, training, eta, 
                kappa, normalizeenergy, normalizestress, descriptors, potential)
    
    return Ae, Af, As, be, bf, bs, nf
end

function ACEdperatomescriptors(x, t, a, b, c, pbc, mldescriptors::Int64)
    return nothing, nothing, nothing, nothing, nothing 
end

function mlperatomdescriptors(x, t, a, b, c, pbc, mldescriptors)
    rcutmax = 0.0
    if mldescriptors.name == "SNAP"         
        for ielem = 1:mldescriptors.ntypes
            rcutmax = max(2.0*mldescriptors.radelem[1+ielem]*mldescriptors.rcutfac,rcutmax)
        end        
        mldescriptors.rcutmax = rcutmax 
        bi, bd, bv, mld, Dd, ai, aj, pairnum = snapdescriptors(x, t, a, b, c, pbc, mldescriptors)
    end
    if mldescriptors.name == "ACE" 
        rcutmax = max(rcutmax, mldescriptors.rcut)
        mld, Dd, ai, aj, pairnum = ACEdperatomescriptors(x, t, a, b, c, pbc, mldescriptors)
    end
    Dd = permutedims(Dd, [2 1 3])     

    dim = size(Dd,1)
    M = size(Dd,3)
    N = size(mld,1)
    mldd = zeros(dim,N,N,M)
    for m = 1:M
        mldd[:,:,:,m] = tallypairforce(Dd[:,:,m], ai, aj, pairnum);
    end
       
    ba = sum(mld,dims=1)
    bb = reshape(sum(mldd,dims=3), (dim,N,M))
    if maximum(bi[:] - ba[:]) > 1e-11
        display(size(mld))
        display(size(Dd))
        display(size(bi))
        display(size(bd))    
        display(maximum(bi[:] - ba[:]))
        error("ba here")
    end
    if maximum(bd[:] - bb[:]) > 1e-11
        display(size(mld))
        display(size(Dd))
        display(size(bi))
        display(size(bd))    
        display(maximum(bd[:] - bb[:]))
        display(bb)
        display(bd)
        error("bb here")
    end

    return mld, mldd, rcutmax  
end

function machinelearningperatomdescriptors(x, t, a, b, c, pbc, mldescriptors)

    rcutmax = 0.0
    if length(mldescriptors) > 0
        mld = 1.0
        mldd = 1.0
        for i = 1:length(mldescriptors)       
            md, mdd, rcut = mlperatomdescriptors(x, t, a, b, c, pbc, mldescriptors[i])                                
            rcutmax = max(rcutmax, rcut)
            if i == 1
                mld = 1.0*md
                mldd = 1.0*mdd
            else
                mld = cat(mld, md, dims=2);     
                mldd = cat(mldd, mdd, dims=4);     
            end
        end
        return mld, mldd, rcutmax
    else
        return nothing, nothing 
    end    
end

function peratomdescriptors(x, q, t, a, b, c, pbc, eta, kappa, emdescriptors, mldescriptors)

    rcutmax = 0.0
    if length(emdescriptors) > 0        
        for i = 1:length(emdescriptors)                
            rcutmax = max(rcutmax, emdescriptors[i].rcut)        
        end
        ed, edd = empiricalperatomdescriptors(x, q, t, a, b, c, pbc, rcutmax, eta, kappa, emdescriptors)
        sz = size(edd)
        if length(sz)==4
            edd = reshape(edd, (sz[1]*sz[2]*sz[3], sz[4]))
        else
            edd = reshape(edd, (sz[1]*sz[2]*sz[3], 1))
        end
    end

    if length(mldescriptors) > 0
        mld, mldd, rcut = machinelearningperatomdescriptors(x, t, a, b, c, pbc, mldescriptors)                                
        sz = size(mldd)
        if length(sz)==4
            mldd = reshape(mldd, (sz[1]*sz[2]*sz[3], sz[4]))
        else
            mldd = reshape(mldd, (sz[1]*sz[2]*sz[3], 1))
        end
        rcutmax = max(rcutmax, rcut)
    end

    if (length(emdescriptors) > 0) & (length(mldescriptors) > 0)        
        d = cat(ed, mld, dims=2);    
        dd = cat(edd, mldd, dims=2);    
        return d, dd, rcutmax 
    elseif length(emdescriptors) > 0
        return ed, edd, rcutmax 
    elseif length(mldescriptors) > 0
        return mld, mldd, rcutmax
    else
        error("empirical descriptors and ML descriptors are empty")
    end    
end

function peratomdescriptors(config, indices, eta, kappa, nspecies, descriptors, potential=nothing)

# getdescriptors
emdescriptors, mldescriptors = getemldescriptors(descriptors)

# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

ne = zeros(Int32, length(indices))
nf = zeros(Int32, length(indices))
be = []
bf = []
Ae = []
Af = []
AA = []
PP = []
for i = 1:length(indices)
    ci = indices[i]
    if ci > nconfigs
        error("Indices are greater than nconfigs")
    end

    x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
    t = config.t[:,(natom[ci]+1):natom[ci+1]] # atom types of configuration ci
    if length(config.q)>0
        q = config.q[:,(natom[ci]+1):natom[ci+1]]
    else
        q = []
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

    d, dd, rcutmax = peratomdescriptors(x, q, t, a[:], b[:], c[:], pbc[:], eta, kappa, emdescriptors, mldescriptors)        
    
    # Adjacency matrix 
    y, alist, neighlist, neighnum = fullneighborlist(x, a[:], b[:], c[:], pbc[:], rcutmax)    
    natoms = length(neighnum)-1
    A = zeros(natoms,natoms)   
    P = zeros(natoms,nspecies)     
    for ii = 1:natoms
        n1 = neighnum[ii];
        n2 = neighnum[ii+1];    
    
        itype = t[ii];
        gg = neighlist[(n1+1):n2]; # a list of neighbors around atom i            
        jj = alist[gg];    # a list of neighbors around atom i            
        #jtype = t[j]; # types of neighbors around atom i
        A[jj,ii] .= 1.0
        A[ii,jj] .= 1.0
        P[ii,itype] = 1.0
    end    
    # display(A)
    # display(P)
    # error("here")

    if (potential !== nothing)  
        potential = Preprocessing.deletenothing(potential)
        if (length(potential) > 0)             
            rcutmax = 0.0
            for j = 1:length(potential)                
                rcutmax = max(rcutmax, potential[j].rcut)            
            end
            eref, e, fref = empiricalpotential(x, q, t, a[:], b[:], c[:], pbc[:], rcutmax, eta, kappa, potential)        
        else
            eref = 0.0; fref = [0.0];
        end
    else
        eref = 0.0; fref = [0.0]; 
    end

    ne[i] = size(d,1)
    nf[i] = size(dd,1)
    if i == 1
        Ae = 1.0*d
        Af = 1.0*dd
        AA = 1.0*A[:]
        PP = 1.0*P[:]
    else
        Ae = cat(Ae, d, dims=1)
        Af = cat(Af, dd, dims=1)
        AA = cat(AA, A[:], dims=1)
        PP = cat(PP, P[:], dims=1)
    end

    if length(config.e)>0
        e = config.e[indices[i]] - eref
        be = [be; e]
    end    
    if length(config.f)>0
        f = config.f[:,(natom[ci]+1):natom[ci+1]]
        bf = [bf; f[:] .- fref[:]]    
    end    
end

return Ae, Af, be, bf, AA, PP, ne, nf

end

function peratomdescriptors(data, descriptors, options, potential=nothing)
        
    data = Preprocessing.deletenothing(data)
    descriptors = Preprocessing.deletenothing(descriptors)
    nspecies = length(data[1].atomspecies)    

    eta = options.eta 
    kappa = options.kappa 
    pbc = options.pbc 
    a = options.a
    b = options.b
    c = options.c
    
    if pbc === nothing        
        pbc = [1; 1; 1] # periodic boundary conditions
    end
    
    n = length(data)
    Ae = 1.0; Af = 1.0; be = 1.0; bf = 1.0; A = 1.0; P = 1.0; 
    ne = Array{Any}(nothing, n);
    nf = Array{Any}(nothing, n);
    for i = 1:n
        # training data set
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        training = Array(1:config.nconfigs)      

        Aei, Afi, bei, bfi, Ai, Pi, nei, nfi = peratomdescriptors(config, training, eta, 
                kappa, nspecies, descriptors, potential)
        
        ne[i] = nei     
        nf[i] = nfi 
        if (i == 1) 
            Ae = 1.0*Aei
            Af = 1.0*Afi
            be = 1.0*bei
            bf = 1.0*bfi
            A = 1.0*Ai
            P = 1.0*Pi 
        else         
            Ae, be = catarray(Ae, be, Aei, bei)
            Af, bf = catarray(Af, bf, Afi, bfi)
            A = [A; Ai]
            P = [P; Pi]
        end                        
    end

    return Ae, Af, be, bf, A, P, ne, nf
end


function peratomdescriptors(data::Preprocessing.DataStruct, descriptors, options, potential=nothing)
        
    descriptors = Preprocessing.deletenothing(descriptors)
    nspecies = length(data.atomspecies)    

    eta = options.eta 
    kappa = options.kappa 
    pbc = options.pbc 
    a = options.a
    b = options.b
    c = options.c
    
    if pbc === nothing        
        pbc = [1; 1; 1] # periodic boundary conditions
    end
    
    config, indices = Preprocessing.readconfigdata(data, pbc, a, b, c)                  
    training = Array(1:config.nconfigs)      

    Ae, Af, be, bf, A, P, ne, nf = peratomdescriptors(config, training, eta, 
            kappa, nspecies, descriptors, potential)

    return Ae, Af, be, bf, A, P, ne, nf
end
