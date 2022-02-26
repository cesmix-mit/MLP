#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module Potential

using Revise, LinearAlgebra
using Base: @kwdef

export empiricalpotential, snappotential, addpotential, getpotential, initsna, initsnap, PotentialStruct
export empiricaldescriptors, snapdescriptors, emldescriptors, DescriptorsOptions, DNNdata, GNNdata, SNAPparams 
export boundingbox, domainmaps, atomlist, neighborlist, fullneighborlist, neighsingles, accumarray 
export neighpairlist, neighpairs, neightripletlist, neightriplets, neighquadletlist, neighquadlets 
export tallysingle, tallypair, tallytriplet, tallyquadlet, findatomtype, periodicimages
export initPOD, POD, PODdesc, PODdescriptors 
#export lammpspotential, lammpssnapdescriptors

@kwdef mutable struct DescriptorsOptions
    normalizeenergy::Bool = true   
    normalizestress::Bool = true   
    pbc = nothing
    a = nothing
    b = nothing
    c = nothing
    eta = nothing 
    kappa = nothing 
end

function catarray(Ae, be, Aei, bei)
    if (length(Aei) > 0) & (length(bei) > 0)
        Ae = cat(Ae, Aei, dims=1)
        be = [be; bei]
    end
    return Ae, be
end

function latticevolume(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64})
    vol = dot(a, cross(b, c))    
    return vol
end

function getemldescriptors(descriptors)

    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   

    n = length(descriptors) 
    emdescriptors = Array{Any}(nothing, n)
    mldescriptors = Array{Any}(nothing, n)

    n = 0
    m = 0
    for i = 1:length(descriptors)    
        if typeof(descriptors[i]) == Potential.PotentialStruct
            n = n + 1
            emdescriptors[n] = descriptors[i]                         
        else
            m = m + 1
            mldescriptors[m] = descriptors[i]                         
        end
    end

    if n>0
        if length(descriptors) > n
            for i = length(descriptors):-1:(n+1)    
                deleteat!(emdescriptors, i)            
            end   
        end     
    else
        emdescriptors = []
    end

    if m>0
        if length(descriptors) > m
            for i = length(descriptors):-1:(m+1)    
                deleteat!(mldescriptors, i)            
            end   
        end     
    else
        mldescriptors = []
    end

    return emdescriptors, mldescriptors
end

include("accumarray.jl");
include("atomlist.jl");
include("boundingbox.jl");
include("domainmaps.jl");
include("findatomtype.jl");
include("periodicimages.jl");
include("prebox.jl");
include("potentialstruct.jl");
include("latticecoords.jl");
include("neighborlist.jl");
include("fullneighborlist.jl");
include("neighsingles.jl");
include("neighpairlist.jl");
include("neighpairs.jl");
include("neightripletlist.jl");
include("neightriplets.jl");
include("neighquadletlist.jl");
include("neighquadlets.jl");
include("empiricalpotential.jl");
include("empiricaldescriptors.jl");
include("tally.jl");
include("podcpp.jl");
include("PODdescriptors.jl");
include("snastruct.jl");
include("snapcpp.jl");
#include("lammpspotential.jl");
#include("lammpssnapdescriptors.jl");
include("emlpotential.jl");
include("descriptors.jl");
include("emldescriptors.jl");
include("savedata.jl");

end

# function getpotentials(potentials)

#     for i = length(potentials):-1:1    
#         if potentials[i] === nothing
#             deleteat!(potentials, i)            
#         end
#     end   

#     n = length(potentials) 
#     empotentials = Array{Any}(nothing, n)
#     mlpotentials = nothing;

#     n = 0
#     for i = 1:length(potentials)    
#         if typeof(potentials[i]) == Potential.PotentialStruct
#             n = n + 1
#             empotentials[n] = potentials[i]                         
#         else
#             mlpotentials = potentials[i]                         
#         end
#     end

#     if n>0
#         if length(potentials) > n
#             for i = length(potentials):-1:(n+1)    
#                 deleteat!(empotentials, i)            
#             end   
#         end     
#     else
#         empotentials = []
#     end

#     return empotentials, mlpotentials
# end



# mutable struct Atoms{T <: AbstractFloat} <: AbstractAtoms{T}
#     X::Vector{JVec{T}}              # positions
#     P::Vector{JVec{T}}              # momenta (or velocities?)
#     M::Vector{T}                    # masses
#     Z::Vector{AtomicNumber}         # atomic numbers
#     cell::JMat{T}                   # cell
#     pbc::JVec{Bool}                 # boundary condition
#     calc::Union{Nothing, AbstractCalculator}
#     dofmgr::DofManager
#     data::Dict{Any,JData{T}}
#  end
