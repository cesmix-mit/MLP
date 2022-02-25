cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "installation/setup.jl");

push!(LOAD_PATH, MDPpath * "src/Julia/ACEpot");
using ACEpot
using DelimitedFiles

# path to the database 
datapath = "../../data/Si/"
dataformat = "jsonpymatgen"
fileextension = "json"
atomspecies = ["Si"];

# weights for energies, forces, and stresses in the linear fit
weightinner=[100  1  1.00E-8]

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 1.0E-8]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = false 

# training data 
traindata[1] = adddata(datapath * "training", dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)

# testing data 
testdata[1] = adddata(datapath * "test", dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
                rotationmatrix, transposelattice)

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# cut-off radius
rcut = 5.0

# inner radius 
rin = 0.9

# bounds for cut-off radius to be optimized
rcutrange = [3.9, 5.0]

# bounds for inner radius to be optimized
rinrange = [0.8, 1.4]

# define range for nonlinear parameters to be optimized
etarange =  hcat(rcutrange, rinrange)'

# number of interpolation points for each nonlinear parameters 
N = [7,7]
etaspace = [etarange N]

# optimization parameters
eta = [rcut, 1.0]; 
kappa = [0];  
optim = setoptim(lossfunc, method, eta, kappa, weightouter, etaspace)

# # ZBL reference potential
# potentials[1] = addpotential("pair", "nonbonded", [1, 1], "zbl", [4.0, 4.8, 73], 4.8);
# include(potentials[1].potentialfunction * ".jl");  # include the potential file

r0 = ACEpot.rnn(:Si)
rin = 0.6*r0 
rcut = 5.0;

# single descriptors to catpure energy of isolated pure elements
descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = ACEpot.ACEparams(species = [:Si], nbody=4, pdegree=12, r0=r0, rcut=rcut, rin=rin, wL=1.5, csp=1.5)

# linear fit to compute SNAP coefficients
coeff, ce, cf, cef, cefs = linearfit(traindata, descriptors, potentials, Doptions, optim)

# compute unweighted MAE, RMSE, RSQ errors 
energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

printerrors(["traing"], energyerrors, "Energy Errors")
printerrors(["traing"], forceerrors, "Force Errors")

# compute unweighted MAE, RMSE, RSQ errors 
energytesterrors, forcetesterrors, stresstesterrors = validate(testdata, descriptors, potentials, Doptions, optim, coeff)    

printerrors(["test"], energytesterrors, "Energy Errors")
printerrors(["test"], forcetesterrors, "Force Errors")

