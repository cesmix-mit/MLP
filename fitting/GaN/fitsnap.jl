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

wj = [1.0, 0.9590493408]
radelem = [0.5, 0.417932464]
bzeroflag = 1
chemflag = 0

# cut-off radius
eta = zeros(5,1)
for j = 1:5
    tm = readdlm("results/optsnapeta" * string(j-1) *  ".txt")
    eta[j] = tm[1]
end

twojmaxlist = [2, 4, 6, 8, 10]
for j = 0:4
    rcut = eta[j+1]
    twojmx = twojmaxlist[j+1]
    display([j twojmx rcut])

    # SNAP descriptors
    descriptors[1] = SNAPparams(species = [:Ga,:N], twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363, 
        elemradius = radelem, elemweight = wj, bzeroflag = bzeroflag, chemflag=chemflag)

    coeff = linearfit(traindata, descriptors, Doptions)
    energyerrors, forceerrors = Potential.snaperroranalysis(traindata, descriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitsnapcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitsnaptrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
end

