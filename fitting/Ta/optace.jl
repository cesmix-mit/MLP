cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

push!(LOAD_PATH, MLPpath * "src/ACEpot");
using ACEpot
using DelimitedFiles

# path to the database 
datapath = "../../data/Ta/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["Ta"];

# folders in the datapath 
folders = ["Displaced_A15", "Displaced_BCC", "Displaced_FCC", "Elastic_BCC",
            "Elastic_FCC", "GSF_110", "GSF_112", "Liquid", "Surface",
            "Volume_A15", "Volume_BCC", "Volume_FCC"];

# weights for energies, forces, and stresses in the linear fit
weightinner=[100            1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09]

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 1.0E-9]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = true 

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

# cut-off radius
rcut = 5.0

# inner radius 
rin = 0.9

# bounds for cut-off radius to be optimized
rcutrange = [4.25, 5.0]

# bounds for inner radius to be optimized
rinrange = [0.9, 1.5]

# define range for nonlinear parameters to be optimized
etarange =  hcat(rcutrange, rinrange)'

# number of interpolation points for each nonlinear parameters 
N = [7,7]
etaspace = [etarange N]

# optimization parameters
eta = [rcut, 1.0]; 
kappa = [0];  
optim = setoptim(lossfunc, method, eta, kappa, weightouter, etaspace)

for j = 3:3

    local acedescriptors = Array{Any}(nothing, 3)

    # single acedescriptors to catpure energy of isolated pure elements
    acedescriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
    
    # ACE acedescriptors
    if j == 0
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=3, pdegree=6, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=8, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=10, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=2, pdegree=12, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=13, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end
    
    # optimize the POD potential 
    opteta, coeff, fmin, iter, polycoeff = optimize(traindata, traindata, acedescriptors, potentials, Doptions, optim)

    # compute unweighted MAE, RMSE, RSQ errors 
    energyerrors, forceerrors, stresserrors = validate(traindata, acedescriptors, potentials, Doptions, optim, coeff)    

    printerrors(folders, energyerrors, "Energy Errors")
    printerrors(folders, forceerrors, "Force Errors")

    # # compute unweighted MAE, RMSE, RSQ errors 
    # energytesterrors, forcetesterrors, stresstesterrors = validate(testdata, acedescriptors, potentials, Doptions, optim, coeff)    

    # printerrors(["test"], energytesterrors, "Energy Errors")
    # printerrors(["test"], forcetesterrors, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optacecoeff" * string(j) *  ".txt", coeff)
    writedlm("results/optacepolycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    #writedlm("results/optacetesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optaceeta" * string(j) *  ".txt", opteta)
end


