module ACEpot 

using ACE1 

using Base: @kwdef

# function rpi_basis(; species = :X, N = 3,
#     # transform parameters
#     r0 = 2.5,
#     trans = PolyTransform(2, r0),
#     # degree parameters
#     D = SparsePSHDegree(),
#     maxdeg = 8,
#     # radial basis parameters
#     rcut = 5.0,
#     rin = 0.5 * r0,
#     pcut = 2,
#     pin = 0,
#     constants = false,
#     rbasis = transformed_jacobi(get_maxn(D, maxdeg, species), trans, rcut, rin;
#                                 pcut=pcut, pin=pin),
#     # one-particle basis
#     Basis1p = BasicPSH1pBasis,
#     basis1p = Basis1p(rbasis; species = species, D = D) )

@kwdef mutable struct ACEparams
    name::String = "ACE"
    species::Vector{Symbol} = [:X]
    nbody::Int64 = 2        
    pdegree::Int64 = 8         
    r0::Float64=2.5  
    rcut::Float64=5.0 
    rin::Float64=1.0
    pcut::Int64=2 
    pin::Int64=0  
    wL::Float64 = 1.5
    csp::Float64 = 1.0    
end

function ACEatoms(x, v, m, z, box, pbc)
    natoms = size(x,2)
    if length(m) == 1
        m = m*ones(natoms)
    end
    if length(z) == 1
        z = z*ones(natoms)
    end
    atom = ACE1.Atoms(X = x, P = v, M = m, Z = z, cell = box, pbc = pbc)    
    return atom
end

function ACEatoms(x, v, m, t, a, b, c, pbc, Z)
    # natoms = size(x,2)
    # if length(m) == 1
    #     m = m*ones(natoms)
    # end
    z = Int32.(0*t) 
    for i = 1:length(Z)
        ind = findall(t .== i)
        z[ind] .= Z[i]
    end
    # Transpose since ACE uses rows as lattice vectors 
    box = [a[:] b[:] c[:]]'                    
    at = ACE1.Atoms(X = x, P = v, M = m, Z = Int32.(z[:]), cell = box, pbc = Bool.(pbc))    
    # display(box)
    # display(at.cell)
    # error("here")
    return at 
end

function ACEbasis(params::ACEparams)
    return ACE1.rpi_basis(;
    species = params.species,
    N       = params.nbody - 1,
    maxdeg  = params.pdegree,
    D       = ACE1.SparsePSHDegree(; wL = params.wL, csp = params.csp),
    r0      = params.r0, 
    rin     = params.rin,
    rcut    = params.rcut,
    pin     = params.pin,
    pcut    = params.pcut,
    )
end

function ACEdescriptors(atoms::Atoms, params::ACEparams)
    basis = ACEbasis(params)
    d = ACE1.energy(basis, atoms)
    dd = ACE1.forces(basis, atoms)
    dv = ACE1.virial(basis, atoms)
    d = d[:]
    M = length(d)
    N = length(dd[1])
    dim = size(dv[1],1)
    df = zeros(dim,N,M)
    ds = zeros(6,M)        
    for i = 1:M
        for j = 1:N
            df[1,j,i] = dd[i][j][1]
            df[2,j,i] = dd[i][j][2]
            df[3,j,i] = dd[i][j][3]
        end
        ds[1,i] = dv[i][1,1]
        ds[2,i] = dv[i][2,2]
        ds[3,i] = dv[i][3,3]
        ds[4,i] = dv[i][3,2]
        ds[5,i] = dv[i][3,1]
        ds[6,i] = dv[i][2,1]
    end                
    return d, df, ds
end

function ACEdescriptors(x, t, a, b, c, pbc, params::ACEparams)       
    natoms = size(x,2)        
    m = ones(natoms)       
    if length(t) == 1
        t = t*ones(natoms)
    end
    Z = zeros(Int32,length(params.species))  
    for i = 1:length(Z)
        Z[i] = AtomicNumber(params.species[i]).z
    end
    atoms = ACEatoms(x, x, m, t, a, b, c, pbc, Z);
    return ACEdescriptors(atoms, params)
end

function ACEdescriptors(x, t, a, b, c, pbc, mldescriptors::Int64)
    return nothing, nothing, nothing
end

function ACEperatomdescriptors(at::Atoms, params::ACEparams)
    shipB = ACEbasis(params)
    D = zeros(fltype(shipB), length(shipB), length(at))
    B = ACE1.alloc_B(shipB)
    nlist = ACE1.neighbourlist(at, ACE1.cutoff(shipB); storelist=false)
    maxnR = ACE1.maxneigs(nlist)
    tmp = ACE1.alloc_temp(shipB, maxnR)
    tmpRZ = (R = zeros(JVec{T}, maxnR), Z = zeros(AtomicNumber, maxnR))

    N =  length(at)
    pairnum = zeros(Int32, N)
    ai = zeros(Int32, N*maxnR)
    aj = zeros(Int32, N*maxnR)
    k = 1
    for i = 1:length(at)
        j, R, Z = ACE1.neigsz!(tmpRZ, nlist, at, i)
        fill!(B, 0)
        ACE1.evaluate!(B, tmp, shipB, R, Z, at.Z[i])
        D[:,i] = B[:]
        pairnum[i] = length(R)
        for a = 1:length(R)
            ai[k] = i
            aj[k] = j[a]
            k = k + 1
        end
    end
    pairnum = [0; cumsum(pairnum)];

    # allocate space accordingly
    Dd = zeros(JVec{T}, length(shipB), maxnR, length(at))
    B = ACE1.alloc_B(shipB, maxnR)
    dB = ACE1.alloc_dB(shipB, maxnR)
    tmp = ACE1.alloc_temp_d(shipB, maxnR)
    tmpRZ = (R = zeros(JVec{T}, maxR), Z = zeros(AtomicNumber, maxnR))
    Dd = ace_forces_inner!(shipB, at, nlist, Dd, B, dB, tmp, tmpRZ)

    return D, Dd, ai, aj, pairnum  
end

function ace_forces_inner!(shipB, at::AbstractAtoms{T},
                        nlist, F, B, dB, tmp, tmpRZ) where {T}
    # assemble site gradients and write into F        
    for i = 1:length(at)
        j, R, Z = ACE1.neigsz!(tmpRZ, nlist, at, i)
        fill!(dB, zero(JVec{T}))
        fill!(B, 0)
        ACE1.evaluate_d!(dB, tmp, shipB, R, Z, at.Z[i])
        F[:,:,1:length(R),i] = dB 
    end
    return F 
end

function ACEperatomdescriptors(x, t, a, b, c, pbc, params::ACEparams)                
    natoms = size(x,2)        
    m = ones(natoms)       
    if length(t) == 1
        t = t*ones(natoms)
    end
    atoms = ACEatoms(x, x, m, t, a, b, c, pbc);
    return ACEperatomdescriptors(atoms, params)
end

end


# function energy(shipB::IPBasis, at::AbstractAtoms{T}) where {T}
#     E = zeros(fltype(shipB), length(shipB))
#     B = alloc_B(shipB)
#     nlist = neighbourlist(at, cutoff(shipB); storelist=false)
#     maxnR = maxneigs(nlist)
#     tmp = alloc_temp(shipB, maxnR)
#     tmpRZ = (R = zeros(JVec{T}, maxnR), Z = zeros(AtomicNumber, maxnR))
#     for i = 1:length(at)
#        j, R, Z = neigsz!(tmpRZ, nlist, at, i)
#        fill!(B, 0)
#        evaluate!(B, tmp, shipB, R, Z, at.Z[i])
#        E[:] .+= B[:]
#     end
#     return E
#  end
 
 
#  function forces(shipB::IPBasis, at::AbstractAtoms{T}) where {T}
#     # precompute the neighbourlist to count the number of neighbours
#     nlist = neighbourlist(at, cutoff(shipB); storelist=false)
#     maxR = maxneigs(nlist)
#     # allocate space accordingly
#     F = zeros(JVec{T}, length(at), length(shipB))
#     B = alloc_B(shipB, maxR)
#     dB = alloc_dB(shipB, maxR)
#     tmp = alloc_temp_d(shipB, maxR)
#     tmpRZ = (R = zeros(JVec{T}, maxR), Z = zeros(AtomicNumber, maxR))
#     return forces_inner!(shipB, at, nlist, F, B, dB, tmp, tmpRZ)
#  end
 
#  # this is a little hack to remove a type instability. It probably makes no
#  # difference in practise...
#  function forces_inner!(shipB::IPBasis, at::AbstractAtoms{T},
#                         nlist, F, B, dB, tmp, tmpRZ) where {T}
#     # assemble site gradients and write into F
#     for i = 1:length(at)
#        j, R, Z = neigsz!(tmpRZ, nlist, at, i)
#        fill!(dB, zero(JVec{T}))
#        fill!(B, 0)
#        evaluate_d!(dB, tmp, shipB, R, Z, at.Z[i])
#        for a = 1:length(R)
#           F[j[a], :] .-= dB[:, a]
#           F[i, :] .+= dB[:, a]
#        end
#     end
#     return [ F[:, iB] for iB = 1:length(shipB) ]
#  end
 




#  function ACEperatom(shipB::IPBasis, at::ACEatoms{T}) where {T}
#     # precompute the neighbourlist to count the number of neighbours
#     nlist = neighbourlist(at, cutoff(shipB); storelist=false)
#     maxR = maxneigs(nlist)
#     # allocate space accordingly
#     F = zeros(JVec{T}, length(shipB), maxR, length(at))
#     B = alloc_B(shipB, maxR)
#     dB = alloc_dB(shipB, maxR)
#     tmp = alloc_temp_d(shipB, maxR)
#     tmpRZ = (R = zeros(JVec{T}, maxR), Z = zeros(AtomicNumber, maxR))
#     return ace_forces_inner!(shipB, at, nlist, F, B, dB, tmp, tmpRZ)
#  end
 
# function ACEdescriptors(atoms::ACEatoms, params::ACEparams)
#     ACE1.site_energy(ACEbasis(params), atoms, )
# end
#function site_energy(basis::IPBasis, at::AbstractAtoms, i0::Integer)

# basis = rpi_basis(; 
#       species = [:Si],
#       N = 3,                        # correlation order = body-order - 1
#       maxdeg = 12,                  # polynomial degree
#       D = SparsePSHDegree(; wL=1.5, csp=1.0),
#       r0 = r0,                      # estimate for NN distance
#       rin = 0.65*r0, rcut = 5.0,    # domain for radial basis (cf documentation)
#       pin = 0)

# basis = rpi_basis(; 
#       species = ["Si","C"],
#       N = 3,                        # correlation order = body-order - 1
#       maxdeg = 12,                  # polynomial degree
#       D = SparsePSHDegree(; wL=1.5, csp=1.0),
#       r0 = r0,                      # estimate for NN distance
#       rin = 0.65*r0, rcut = 5.0,    # domain for radial basis (cf documentation)
#       pin = 0)


# function rpi_basis(; species = :X, N = 3,
#     # transform parameters
#     r0 = 2.5,
#     trans = PolyTransform(2, r0),
#     # degree parameters
#     D = SparsePSHDegree(),
#     maxdeg = 8,
#     # radial basis parameters
#     rcut = 5.0,
#     rin = 0.5 * r0,
#     pcut = 2,
#     pin = 0,
#     constants = false,
#     rbasis = transformed_jacobi(get_maxn(D, maxdeg, species), trans, rcut, rin;
#                                 pcut=pcut, pin=pin),
#     # one-particle basis
#     Basis1p = BasicPSH1pBasis,
#     basis1p = Basis1p(rbasis; species = species, D = D) )

#  return RPIBasis(basis1p, N, D, maxdeg, constants)
# end

