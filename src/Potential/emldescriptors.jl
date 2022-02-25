using Preprocessing

# function ACEdescriptors(x, t, a, b, c, pbc, mldescriptors::Int64)
#     return nothing, nothing, nothing
# end

# function machinelearningdescriptors(x, t, a, b, c, pbc, mldescriptors)
#     if mldescriptors.name == "SNAP" 
#         rcutmax = 0.0
#         for ielem = 1:mldescriptors.ntypes
#             rcutmax = max(2.0*mldescriptors.radelem[1+ielem]*mldescriptors.rcutfac,rcutmax)
#         end        
#         mldescriptors.rcutmax = rcutmax 
#         mld, mldd, mldv = snapdescriptors(x, t, a, b, c, pbc, mldescriptors)
#         sz = size(mldd)
#         mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
#     end
#     # if mldescriptors.name == "ACE" 
#     #     mld, mldd, mldv = ACEdescriptors(x, t, a, b, c, pbc, mldescriptors)
#     #     sz = size(mldd)
#     #     mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
#     # end
#     return mld, mldd, mldv
# end

function emldescriptors(x, q, t, a, b, c, pbc, eta, kappa, emdescriptors, mldescriptors=nothing)

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
    if mldescriptors !== nothing
        if isa(mldescriptors, Array)
            for i = 1:length(mldescriptors)       
                md, mdd, mdv = machinelearningdescriptors(x, t, a, b, c, pbc, mldescriptors[i])                                
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
        else
            if mldescriptors.name == "SNAP" 
                rcutmax = 0.0
                for ielem = 1:mldescriptors.ntypes
                    rcutmax = max(2.0*mldescriptors.radelem[1+ielem]*mldescriptors.rcutfac,rcutmax)
                end        
                mldescriptors.rcutmax = rcutmax 
                mld, mldd, mldv = snapdescriptors(x, t, a, b, c, pbc, mldescriptors)
                sz = size(mldd)
                mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
            end
            # if mldescriptors.name == "ACE" 
            #     mld, mldd, mldv = ACEdescriptors(x, t, a, b, c, pbc, mldescriptors)
            #     sz = size(mldd)
            #     mldd = reshape(mldd, (sz[1]*sz[2], sz[3]))
            # end
        end
    end

    if (length(emdescriptors) > 0) & (mldescriptors !== nothing)
        d = [ed[:]; mld[:]]
        dd = cat(edd, mldd, dims=2);    
        dv = cat(edv, mldv, dims=2);    
        return d, dd, dv    
    elseif length(emdescriptors) > 0
        return ed, edd, edv
    elseif (mldescriptors !== nothing)
        return mld, mldd, mldv 
    else
        error("empirical descriptors and ML descriptors are empty")
    end    
end

function emldescriptors(config, indices, eta, kappa, normalizeenergy, emdescriptors, 
                mldescriptors=nothing, potential=nothing)

normalizestress = 0
# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

nf = zeros(Int64, length(indices))
be = []
bf = []
bs  = []
Ae = []
Af = []
As = []
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

    d, dd, dv = emldescriptors(x, q, t, a[:], b[:], c[:], pbc[:], eta, kappa, emdescriptors, mldescriptors)        

    if normalizestress==1
        vol = latticevolume(a[:], b[:], c[:])
        dv = (1.6021765e6/vol) * dv
    end

    if potential !== nothing        
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
        eref = 0.0;
        fref = [0.0];
        sref = [0.0];
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

function emldescriptors2(config, indices, eta, kappa, normalizeenergy, normalizestress, 
        descriptors, potential=nothing)

# getdescriptors
emdescriptors, mldescriptors = getdescriptors(descriptors)

# cumalative sum of numbers of atoms 
nconfigs = length(config.natom)
natom = [0; cumsum(config.natom[:])];

nf = zeros(Int64, length(indices))
be = []
bf = []
bs  = []
Ae = []
Af = []
As = []
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

    d, dd, dv = emldescriptors(x, q, t, a[:], b[:], c[:], pbc[:], eta, kappa, emdescriptors, mldescriptors)        
    
    if normalizestress==1
        vol = latticevolume(a[:], b[:], c[:])
        dv = (1.6021765e6/vol) * dv
    end

    if (potential !== nothing)  
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

# function catarray(Ae, be, Aei, bei)
#     if (length(Aei) > 0) & (length(bei) > 0)
#         Ae = cat(Ae, Aei, dims=1)
#         be = [be; bei]
#     end
#     return Ae, be
# end

# using Base: @kwdef

# @kwdef mutable struct DescriptorsOptions
#     normalizeenergy::Bool = true   # use gpu (if cuda available)
#     normalizestress::Bool = true   # use gpu (if cuda available)
#     pbc = nothing
#     a = nothing
#     b = nothing
#     c = nothing
#     eta = nothing 
#     kappa = nothing 
# end

function emldescriptors(data, descriptors, options)
        
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

        Aei, Afi, Asi, bei, bfi, bsi, nfi = emldescriptors2(config, training, eta, kappa, 
                normalizeenergy, normalizestress, descriptors)
        
        # display(config)                
        # display(size(Aei))

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

