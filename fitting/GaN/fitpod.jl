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

# Bessel scaling parameters
gamma = [0.0, 2, 4]

# cut-off radius
eta = zeros(5,1)
for j = 1:5
    tm = readdlm("results/optpodeta" * string(j-1) *  ".txt")
    eta[j] = tm[1]
end

# inner radius 
rin = 0.5

for j = 0:4
    rcut = eta[j+1]
    display([j rcut])

    # POD Descriptors
    if j == 0
        descriptors[1] = POD(nbody=2, species = [:Ga,:N], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:Ga,:N], pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j == 1
        # 3*3 + (3*3+3)*6 = 9 + 12*6 = 81
        descriptors[1] = POD(nbody=2, species = [:Ga,:N], pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:Ga,:N], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==2
        descriptors[1] = POD(nbody=2, species = [:Ga,:N], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:Ga,:N], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==3 
        descriptors[1] = POD(nbody=2, species = [:Ga,:N], pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:Ga,:N], pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==4
        descriptors[1] = POD(nbody=2, species = [:Ga,:N], pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
        descriptors[2] = POD(nbody=3, species = [:Ga,:N], pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
    elseif j==5           
        descriptors[1] = POD(nbody=2, species = [:Ga,:N], pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
        descriptors[2] = POD(nbody=3, species = [:Ga,:N], pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
    end

    coeff = linearfit(traindata, descriptors, Doptions)
    energyerrors, forceerrors = Potential.poderroranalysis(traindata, descriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitpodcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitpodtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
end

