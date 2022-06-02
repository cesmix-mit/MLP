cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/HfB2/"
folders = ["EOS1","EOS2","Deformed1","Deformed2","Deformed3","Distorted1","Distorted2","MD1","MD2"]
#folders = ["EOS1","EOS2","Deformed2","Deformed3","Distorted1","Distorted2"]
#folders = ["EOS1","EOS2","Deformed3"]
folders = ["MD1","MD2"]
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Hf", "B"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 100.0
weightinner[:,2] .= 1.0
# weightinner[2,1] = 1e-14
# weightinner[2,2] = 1e-14

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 0.0]

# randomly selecting the configurations in the database
randomize = false;

# use percentage of the data 
percentage = 40.0;

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

# randomly selecting the configurations in the database
#randomize = true;
# test data 
for i = 1:n
testdata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
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

j = 1;
#for j = 0:3
    display([j rcut])

    # POD Descriptors
    if j == 0
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j == 1
        # 3*3 + (3*3+3)*6 = 9 + 12*6 = 81
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-8)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-8)
    elseif j==2
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==3 
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==4
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
    elseif j==5           
        descriptors[1] = POD(nbody=2, species = [:Hf,:B], pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
        descriptors[2] = POD(nbody=3, species = [:Hf,:B], pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
    end

    coeff = linearfit(traindata, descriptors, Doptions)
    energyerrors, forceerrors = Potential.poderroranalysis(testdata, descriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    display(coeff)
    display(e1)
    display(e2)    
    # printerrors(["train"], e1, "Energy Errors")
    # printerrors(["train"], e2, "Force Errors")

    # energytesterrors, forcetesterrors = Potential.poderroranalysis(testdata, descriptors, Doptions, coeff)

    # e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    # e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    # printerrors(["test"], e1, "Energy Errors")
    # printerrors(["test"], e2, "Force Errors")

    # Preprocessing.mkfolder("results")
    # writedlm("results/fitpodcoeff" * string(j) *  ".txt", coeff)
    # writedlm("results/fitpodtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    # writedlm("results/fitpodtesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
#end

