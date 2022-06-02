cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/HfB2/"
folders = ["training"]
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Hf", "B"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 20.0
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

# test data 
testfolders = ["test"]
testdata[1] = adddata(datapath * testfolders[1], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# bounds for cut-off radius to be optimized
rcutrange = [3.25, 5.0]

# define range for nonlinear parameters to be optimized
etarange =  reshape(rcutrange,(1,2))

# number of interpolation points for each nonlinear parameters 
N = [7]
etaspace = [etarange N]

# cut-off radius
rcut = 5.0

# inner radius 
rin = 0.5

# Bessel scaling parameters
gamma = [0.0, 2, 4]

# optimization parameters
eta = [rcut]; 
kappa = [0];  
optim = setoptim(lossfunc, method, eta, kappa, weightouter, etaspace)

j = 4;
#for j = 1:1
    display(j)
    # POD Descriptors
    if j == 0
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-12)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-12)
    elseif j == 1
        # 3*3 + (3*3+3)*6 = 9 + 12*6 = 81
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-10)
        #descriptors[4] = POD(nbody=4, species = [:Hf,:B], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma, b2b4=16,b3b4=24)
    elseif j==2
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-8)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-8)
    elseif j==3 
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
    elseif j==4
        descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-4)
        descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-4)
    elseif j==5           
        descriptors[1] = POD(nbody=2, species = [:Hf,:B], pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
        descriptors[2] = POD(nbody=3, species = [:Hf,:B], pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=true, projectiontol=1e-6)
    end

    # optimize the pod potential 
    opteta, coeff, fmin, iter, polycoeff, etapts, lossvalues, energyerrors, forceerrors = optimize(traindata, testdata, descriptors, potentials, Doptions, optim)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    energytesterrors, forcetesterrors = Potential.poderroranalysis(testdata, descriptors, Doptions, coeff)

    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(["test"], e1, "Energy Errors")
    printerrors(["test"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optpod23Bcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/optpod23Bpolycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optpod23Btrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/optpod23Btesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optpod23Beta" * string(j) *  ".txt", opteta)
#end

