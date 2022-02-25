#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module Optimization

using Base: @kwdef

export setoptim, printerrors, linearfit, optimize, validate
export mpolyeval, mpolyinterp, gradientdescent, num2string, replacestring

@kwdef mutable struct OptimStruct
    lossfunc::String = "energyforce"
    method::String = "lsq"    
    eta::Vector{Float64} = [0.0]
    kappa::Vector{Int64} = [0]
    weightouter::Vector{Float64} = [1, 1, 1]
    etaspace::Matrix{Float64} = [0.0 1.0]
    kappalist::Matrix{Int64}  = [0 0]  
end

function setoptim(lossfunc=nothing, method=nothing, eta=nothing, kappa=nothing, weightouter=nothing, etaspace=nothing, kappalist=nothing)

    if lossfunc === nothing
        lossfunc = "energyforce"
    end
    if method === nothing
        method = "lsq"
    end
    if eta === nothing
        eta = []
    end
    if kappa === nothing
        kappa = []
    end
    if weightouter === nothing
        weightouter = [1.0, 1.0, 1.0]
    end
    if etaspace === nothing
        etaspace = [0.0 1.0]
    end
    if kappalist === nothing
        kappalist = [0 0]
    end

    optim = OptimStruct(lossfunc, method, eta, kappa, weightouter, etaspace, kappalist)

    return optim
end

function getoptim(optim::OptimStruct)

    lossfunc = optim.lossfunc
    method = optim.method    
    eta = optim.eta
    kappa = optim.kappa
    weightouter = optim.weightouter
    etaspace = optim.etaspace
    kappalist = optim.kappalist 

    return lossfunc, method, eta, kappa, weightouter, etaspace, kappalist
end

include("printerrors.jl");
include("tensorpolyfit.jl");
include("linearfit.jl");

#include("OptimStruct.jl");
#include("tensorpolyopt.jl");
#include("maespace.jl");

end





