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

wj = [1.0, 0.9590493408]
radelem = [0.5, 0.417932464]
bzeroflag = 1
chemflag = 0

twojmaxlist = [2, 4, 6, 8, 10]
for j = 0:4
    twojmx = twojmaxlist[j+1]
    display(j)
 
    # SNAP descriptors
    descriptors[1] = SNAPparams(species = [:Hf,:B], nbody = 1, twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363, 
        elemradius = radelem, elemweight = wj, bzeroflag = bzeroflag, chemflag=chemflag)
    descriptors[2] = SNAPparams(species = [:Hf,:B], twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363, 
        elemradius = radelem, elemweight = wj, bzeroflag = bzeroflag, chemflag=chemflag)
    
    # optimize the pod potential 
    opteta, coeff, fmin, iter, polycoeff, etapts, lossvalues, energyerrors, forceerrors = optimize(traindata, testdata, descriptors, potentials, Doptions, optim)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    energytesterrors, forcetesterrors = Potential.snaperroranalysis(testdata, descriptors, Doptions, coeff)

    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(["test"], e1, "Energy Errors")
    printerrors(["test"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optsnapBcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/optsnapBpolycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optsnapBtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/optsnapBtesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optsnapBeta" * string(j) *  ".txt", opteta)
end


