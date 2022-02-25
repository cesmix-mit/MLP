cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/Mo/"
dataformat = "jsonpymatgen"
fileextension = "json"
atomspecies = ["Mo"];
folders = ["training"];

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
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
                rotationmatrix, transposelattice)
end

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
rcutrange = [4.5, 5.5]

# define range for nonlinear parameters to be optimized
etarange =  reshape(rcutrange,(1,2))

# number of interpolation points for each nonlinear parameters 
N = [7]
etaspace = [etarange N]

# optimization parameters
eta = [rcut]; 
kappa = [0];  
optim = setoptim(lossfunc, method, eta, kappa, weightouter, etaspace)

# single descriptors to catpure energy of isolated pure elements
descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)

twojmaxlist = [2, 4, 6, 8, 10, 12]
for j = 3:3
    twojmx = twojmaxlist[j+1]

    # SNAP descriptors
    descriptors[2] = SNAPparams(species = [:Mo], twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363)

    # optimize the SNAP potential 
    opteta, coeff, fmin, iter, polycoeff = optimize(traindata, traindata, descriptors, potentials, Doptions, optim)

    # compute unweighted MAE, RMSE, RSQ errors 
    energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(["train"], energyerrors, "Energy Errors")
    printerrors(["train"], forceerrors, "Force Errors")

    # compute unweighted MAE, RMSE, RSQ errors 
    energytesterrors, forcetesterrors, stresstesterrors = validate(testdata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(["test"], energytesterrors, "Energy Errors")
    printerrors(["test"], forcetesterrors, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optsnapcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/optsnappolycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optsnaptrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/optsnaptesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optsnapeta" * string(j) *  ".txt", opteta)
end


