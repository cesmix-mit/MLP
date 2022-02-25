cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setup.jl");

# path to the database 
datapath = "../../data/Si/"
dataformat = "jsonpymatgen"
fileextension = "json"
atomspecies = ["Si"];
folders = ["training"];

# weights for energies, forces, and stresses in the linear fit
weightinner=[100            1               1.00E-8]

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

rcut = 5.0
rin = 0.9

descriptors[1] = POD(nbody=1, atomtypes=[1], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)

# descriptors[2] = POD(nbody=2, atomtypes=[1,1], pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut)
# descriptors[3] = POD(nbody=3, atomtypes=[1,1,1], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut)
# descriptors[2] = POD(nbody=2, atomtypes=[1,1], pdegree=[4,10], nbasis = [12], rin = rin, rcut=rcut, gamma0 = [0.0, 2, 4])
# descriptors[3] = POD(nbody=3, atomtypes=[1,1,1], pdegree=[4,8,5], nbasis = [12, 5], rin = rin, rcut=rcut, gamma0 = [0.0, 2, 4])

descriptors[2] = POD(nbody=2, pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = [0.0, 2, 4])
descriptors[3] = POD(nbody=3, pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = [0.0, 2, 4])

# descriptors[2] = initPOD(PODdesc(nbody=2, atomtypes=[1,1], pdegree=[4,6], nbasis = [7], rin = 1.0, rcut=rcut))
# descriptors[3] = initPOD(PODdesc(nbody=3, atomtypes=[1,1,1], pdegree=[4,6,6], nbasis = [7, 6], rin = 1.0, rcut=rcut))

# single descriptors to catpure energy of isolated pure elements
# descriptors[1] = PotentialStruct("single", "nonbonded", [1], "single", [0.0], rcut);
# include(descriptors[1].potentialfunction * ".jl");  # include the descriptors file      
# descriptors[2] = initPOD(PODdesc(nbody=2, atomtypes=[1,1], pdegree=[6,12], nbasis = [13], rin = rin, rcut=rcut))
# descriptors[3] = initPOD(PODdesc(nbody=3, atomtypes=[1,1,1], pdegree=[6,12,15], nbasis = [8,15], rin = rin, rcut=rcut))

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# optimization parameters
optim = setoptim(lossfunc, method)

# linear fit to compute SNAP coefficients
coeff, ce, cf, cef, cefs = linearfit(traindata, descriptors, potentials, Doptions, optim)

print("Coeffficients: "), show(stdout, "text/plain", coeff)

# compute unweighted MAE, RMSE, RSQ errors 
energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

printerrors(folders, energyerrors, "Energy Errors")
printerrors(folders, forceerrors, "Force Errors")
printerrors(folders, stresserrors, "Stress Errors")

err = [energyerrors forceerrors stresserrors]
using DelimitedFiles
Preprocessing.mkfolder("results")
writedlm("results/fitpodcoeff.txt", coeff)
writedlm("results/fitpoderror.txt", err)

