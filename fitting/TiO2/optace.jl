cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACEpot

# path to the database 
datapath = "../../data/"
folders = ["TiO2"]
dataformat = "ann"
fileextension = "xsf"
atomspecies = ["Ti", "O"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 100.0
weightinner[:,2] .= 1.0

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 0.0]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 50.0;

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


# randomly selecting the configurations in the database
randomize = true;            
# test data 
testdata[1] = adddata(datapath * folders[1], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# bounds for cut-off radius to be optimized
rcutrange = [3.5, 5.5]

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

for j = 2:4
    display(j)

    local acedescriptors = Array{Any}(nothing, 3)
    # ACE descriptors
    if j == 0
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ti,:O], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ti,:O], nbody=3, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ti,:O], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ti,:O], nbody=4, pdegree=8, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ti,:O], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ti,:O], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ti,:O], nbody=4, pdegree=10, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ti,:O], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ti,:O], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ti,:O], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ti,:O], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ti,:O], nbody=2, pdegree=12, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ti,:O], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        acedescriptors[1] = ACEpot.ACEparams(species = [:Ti,:O], nbody=1, pdegree=6, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ti,:O], nbody=4, pdegree=13, r0=ACEpot.rnn(:Ti), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end

    # optimize the pod potential 
    opteta, coeff, fmin, iter, polycoeff, etapts, lossvalues, energyerrors, forceerrors = optimize(traindata, traindata, acedescriptors, potentials, Doptions, optim)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    energytesterrors, forcetesterrors = Potential.aceerroranalysis(testdata, acedescriptors, Doptions, coeff)

    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(["test"], e1, "Energy Errors")
    printerrors(["test"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optacecoeff" * string(j) *  ".txt", coeff)
    writedlm("results/optacepolycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/optacetesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optaceeta" * string(j) *  ".txt", opteta)
end

