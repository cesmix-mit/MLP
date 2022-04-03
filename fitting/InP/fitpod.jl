cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/InP/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["In", "P"];

# folders in the datapath 
folders = ["aa", "aIn", "aP", "Bulk",
            "EOS", "iIn", "iP", "s_aa", "s_aIn",
            "s_aP", "Shear", "s_iIn", "s_iP", "Strain", 
            "s_vIn","s_vP", "s_vv", "vP", "vv"];
            
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
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
                rotationmatrix, transposelattice)
end

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# optimization parameters
optim = setoptim(lossfunc, method)

# Bessel scaling parameters
gamma = [0.0, 2, 4]
for j = 4:4

    display(j)

    # POD Descriptors
    if j == 0
        rcut = 4.55040786223262  
        rin = 0.3
        descriptors[1] = POD(nbody=2, species = [:In,:P], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
        descriptors[2] = POD(nbody=3, species = [:In,:P], pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
    elseif j == 1
        rcut = 4.2500
        rin = 0.3098365471
        descriptors[1] = POD(nbody=2, species = [:In,:P], pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
        descriptors[2] = POD(nbody=3, species = [:In,:P], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
    elseif j==2
        rcut = 4.1598216734
        rin = 0.359273786
        descriptors[1] = POD(nbody=2, species = [:In,:P], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
        descriptors[2] = POD(nbody=3, species = [:In,:P], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
    elseif j==3 
        rcut = 4.09937356
        rin = 0.397370934
        descriptors[1] = POD(nbody=2, species = [:In,:P], pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
        descriptors[2] = POD(nbody=3, species = [:In,:P], pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
    elseif j==4
        rcut = 3.895327
        rin = 0.4649875
        descriptors[1] = POD(nbody=2, species = [:In,:P], pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
        descriptors[2] = POD(nbody=3, species = [:In,:P], pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma, hybdrid23=true)
    end

    Potential.podprojection(traindata, descriptors, Doptions.pbc)

    # linear fit to compute SNAP coefficients
    coeff, ce, cf, cef, cefs, emae, fmae, smae,~ = linearfit(traindata, descriptors, potentials, Doptions, optim)

    #print("POD Coeffficients: "), show(stdout, "text/plain", coeff)

    # # compute unweighted MAE, RMSE, RSQ erroxrs 
    energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(folders, energyerrors, "Energy Errors")
    printerrors(folders, forceerrors, "Force Errors")

    err = [energyerrors forceerrors]
    display(err)

    using DelimitedFiles
    Preprocessing.mkfolder("results")
    writedlm("results/fitpod23coeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitpod23error" * string(j) *  ".txt", err)
end

