#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module DML 

using Base: @kwdef

@kwdef mutable struct ModelOptions
    η::Array{Float64,1} = [1e-4, 1e-4, 1e-4]# learning rate
    β::Array{Float64,1} = [0.1, 0.9]  # loss function parameters 
    featurenorm::Array{Float64,1} = [1] # normalization of input features 
    batchsize::Int32 = 32       # batch size
    epochs::Int32 = 100         # number of epochs
    regparam::Float64 = 0.0 
    normalizeenergy::Bool = false  
    use_cuda::Bool = false      # use gpu (if cuda available)

    filename::String = "trainbatch"   
    modeltype::String = "dnn"   # machine learning model 
    num_inputs::Int32 = 64      # number of inputs in the convolution layers  
    num_layers::Int32 = 2       # number of convolution layers     
    f::Function                 # activation function for convolution layers     
    df::Function                # derivative of activation function for convolution layers     
    lossfunc::Function = lossmse # default loss function 
    initW::Function = zeros     # initialization function for weight matrices  
    initb::Function = zeros     # initialization function for bias vectors 
    initc::Function = ones     # initialization function for output weights     
    initw::Function = ones     # initialization function 
    inits::Function = ones     # initialization function 
end

include("blas.jl");
include("loaddata.jl");
include("lossfunctions.jl");
include("DNNmodel.jl");
include("DiffNNmodel.jl");
include("GNNmodel.jl");
include("train.jl");

end


