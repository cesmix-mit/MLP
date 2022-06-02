cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACEpot

# path to the database 
datapath = "../../data/HfB2/training/"
folders = ["n3","n6","n24","n48","n81"]
#folders = ["n24","n81","n48","n6","n3"]
#folders = ["EOS1","EOS2","Deformed3"]
#folders = ["MD1","MD2"]
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Hf", "B"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 20.0
weightinner[:,2] .= 1.0
# weightinner[2,1] = 1e-14
# weightinner[2,2] = 1e-14

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 0.0]

# randomly selecting the configurations in the database
randomize = false;

# use percentage of the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors 
transposelattice = false 

# training data 
for i = 1:n
traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
            rotationmatrix, transposelattice)
end

# for i = 1:length(folders)
# testdata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
#             percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
#             rotationmatrix, transposelattice)
# end

testpath = "../../data/HfB2/test/"
testfolders = ["g3","g6"]
for i = 1:length(testfolders)
testdata[i] = adddata(testpath * testfolders[i], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)
end

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# Bessel scaling parameters
gamma = [0.0, 2, 4]

# cut-off radius
rcut = 5.0

# inner radius 
rin = 0.5

for j = 0:4
    display([j rin rcut])

    local acedescriptors = Array{Any}(nothing, 3)
    # ACE descriptors
    if j == 0
        acedescriptors[1] = ACEpot.ACEparams(species = [:Hf,:B], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Hf,:B], nbody=3, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        acedescriptors[1] = ACEpot.ACEparams(species = [:Hf,:B], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Hf,:B], nbody=4, pdegree=8, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2
        acedescriptors[1] = ACEpot.ACEparams(species = [:Hf,:B], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Hf,:B], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Hf,:B], nbody=4, pdegree=10, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        acedescriptors[1] = ACEpot.ACEparams(species = [:Hf,:B], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Hf,:B], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Hf,:B], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        acedescriptors[1] = ACEpot.ACEparams(species = [:Hf,:B], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Hf,:B], nbody=2, pdegree=12, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Hf,:B], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        acedescriptors[1] = ACEpot.ACEparams(species = [:Hf,:B], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Hf,:B], nbody=4, pdegree=13, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end

    coeff = readdlm("results/fitacecoeff" * string(j) *  ".txt")            
    #coeff = linearfit(traindata, acedescriptors, Doptions)

    # energyerrors, forceerrors = Potential.aceerroranalysis(traindata, acedescriptors, Doptions, coeff)

    # e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    # e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    # printerrors(folders, e1, "Energy Errors")
    # printerrors(folders, e2, "Force Errors")

    # Preprocessing.mkfolder("results")
    # writedlm("results/fitacecoeff" * string(j) *  ".txt", coeff)
    # writedlm("results/fitacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    

    # energytesterrors, forcetesterrors = Potential.poderroranalysis(testdata, descriptors, Doptions, coeff)
    # e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    # e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    # printerrors(["test"], e1, "Energy Errors")
    # printerrors(["test"], e2, "Force Errors")

    # Preprocessing.mkfolder("results")
    # writedlm("results/fitpodcoeff" * string(j) *  ".txt", coeff)
    # writedlm("results/fitpodtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    # writedlm("results/fitpodtesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    

    energytesterrors, forcetesterrors, stresserrors, eerr, ferr, ferm, energies = Potential.aceerroranalysis(testdata, acedescriptors, Doptions, coeff)
    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(testfolders, e1, "Energy Errors")
    printerrors(testfolders, e2, "Force Errors")
    writedlm("results/fitacetesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/fitaceestenergies" * string(j) *  ".txt", [energies[1]*3 energies[2]*6])    
end

