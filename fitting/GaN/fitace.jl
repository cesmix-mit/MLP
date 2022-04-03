cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/"
folders = ["GaN"]
dataformat = "extxyz"
fileextension = "exyz"
atomspecies = ["Ga", "N"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 100.0
weightinner[:,2] .= 1.0

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 0.0]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors 
transposelattice = false 

# training data 
traindata[1] = adddata(datapath * folders[1], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# cut-off radius
eta = zeros(5,1)
for j = 1:5
    tm = readdlm("results/optaceeta" * string(j-1) *  ".txt")
    eta[j] = tm[1]
end

# inner radius 
rin = 0.5

for j = 0:4
    rcut = eta[j+1]
    display([j rcut])

    local acedescriptors = Array{Any}(nothing, 3)
    # ACE descriptors
    if j == 0
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=3, pdegree=6, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=4, pdegree=8, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=2, pdegree=3, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        #acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=3, pdegree=4, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ga,:N], nbody=4, pdegree=10, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=2, pdegree=3, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        #acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=3, pdegree=4, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ga,:N], nbody=4, pdegree=12, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=2, pdegree=12, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ga,:N], nbody=4, pdegree=12, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ga,:N], nbody=4, pdegree=13, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end

    coeff = linearfit(traindata, acedescriptors, Doptions)
    energyerrors, forceerrors = Potential.aceerroranalysis(traindata, acedescriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitacecoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
end

