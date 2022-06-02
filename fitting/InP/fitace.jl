cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

push!(LOAD_PATH, MLPpath * "src/ACEpot");
using ACEpot
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

for j = 2:2
    display(j)

    # ACE descriptors
    if j == 0
        descriptors[1] = ACEpot.ACEparams(species = [:In,:P], nbody=3, pdegree=6, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        descriptors[1] = ACEpot.ACEparams(species = [:In,:P], nbody=4, pdegree=8, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2        
        # cut-off radius
        rcut = 4.90866  
        # inner radius 
        rin = 1.23987
        descriptors[1] = ACEpot.ACEparams(species = [:In,:P], nbody=2, pdegree=3, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        descriptors[2] = ACEpot.ACEparams(species = [:In,:P], nbody=4, pdegree=10, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        # 4.9086554390 1.408655439 0.0008903618 0.01622093
        descriptors[1] = ACEpot.ACEparams(species = [:In,:P], nbody=2, pdegree=3, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        descriptors[2] = ACEpot.ACEparams(species = [:In,:P], nbody=4, pdegree=12, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        # 4.9086554390 1.408655439 0.0006962715 0.01637712
        descriptors[1] = ACEpot.ACEparams(species = [:In,:P], nbody=2, pdegree=12, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        descriptors[2] = ACEpot.ACEparams(species = [:In,:P], nbody=4, pdegree=12, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        descriptors[1] = ACEpot.ACEparams(species = [:In,:P], nbody=4, pdegree=13, r0=ACEpot.rnn(:Be), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end
        
    coeff, ce, cf, cef, cefs, emae, fmae, smae,~ = linearfit(traindata, descriptors, potentials, Doptions, optim)
    
    # compute unweighted MAE, RMSE, RSQ errors 
    energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(folders, energyerrors, "Energy Errors")
    printerrors(folders, forceerrors, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitacecoeff" * string(j) *  ".txt", coeff)    
    writedlm("results/fitacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
end


