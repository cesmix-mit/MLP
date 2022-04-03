# Add MDP to Julia search path
push!(LOAD_PATH, MLPpath * "src/Preprocessing");
push!(LOAD_PATH, MLPpath * "src/Potential");
push!(LOAD_PATH, MLPpath * "src/Optimization");
push!(LOAD_PATH, MLPpath * "src/DML");

# Set Julia's PATH enviroment variable so that MDP can call external programs
ENV["PATH"] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";

using Revise, Preprocessing, Potential, Optimization, DML  

# load Snap library 
Potential.loadsnap(MLPpath)

# load POD library 
Potential.loadpod(MLPpath)

# initialize arrays with nothing 
traindata = Array{Any}(nothing, 100)
testdata = Array{Any}(nothing, 100)
validdata = Array{Any}(nothing, 100)
descriptors = Array{Any}(nothing, 100)
potentials  = Array{Any}(nothing, 100)
