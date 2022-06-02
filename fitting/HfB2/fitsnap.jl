cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/HfB2/training/"
folders = ["n3","n6","n24","n48","n81"]
#folders = ["n24","n81","n48","n6","n3"]
#folders = ["EOS1","EOS2","Deformed3"]
#folders = ["MD1","MD2"]
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Hf", "B"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 20.0
weightinner[:,2] .= 1.0
# weightinner[2,1] = 1e-14
# weightinner[2,2] = 1e-14

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 0.0]

# randomly selecting the configurations in the database
randomize = false;

# use percentage of the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors 
transposelattice = false 

# training data 
for i = 1:n
traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
            rotationmatrix, transposelattice)
end

# for i = 1:length(folders)
# testdata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
#             percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
#             rotationmatrix, transposelattice)
# end

testpath = "../../data/HfB2/test/"
testfolders = ["g3","g6"]
for i = 1:length(testfolders)
testdata[i] = adddata(testpath * testfolders[i], dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)
end

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# Bessel scaling parameters
gamma = [0.0, 2, 4]

# cut-off radius
rcut = 5.0
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

    coeff = readdlm("results/fitsnapcoeff" * string(j) *  ".txt")            
    #coeff = linearfit(traindata, descriptors, Doptions)

    energyerrors, forceerrors = Potential.snaperroranalysis(traindata, descriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(folders, e1, "Energy Errors")
    printerrors(folders, e2, "Force Errors")

    # Preprocessing.mkfolder("results")
    # writedlm("results/fitsnapcoeff" * string(j) *  ".txt", coeff)
    # writedlm("results/fitsnaptrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    

    # energytesterrors, forcetesterrors = Potential.poderroranalysis(testdata, descriptors, Doptions, coeff)
    # e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    # e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    # printerrors(["test"], e1, "Energy Errors")
    # printerrors(["test"], e2, "Force Errors")

    # Preprocessing.mkfolder("results")
    # writedlm("results/fitpodcoeff" * string(j) *  ".txt", coeff)
    # writedlm("results/fitpodtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    # writedlm("results/fitpodtesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    

    energytesterrors, forcetesterrors, stresserrors, eerr, ferr, ferm, energies = Potential.snaperroranalysis(testdata, descriptors, Doptions, coeff)
    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(testfolders, e1, "Energy Errors")
    printerrors(testfolders, e2, "Force Errors")
    writedlm("results/fitsnaptesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/fitsnapestenergies" * string(j) *  ".txt", [energies[1]*3 energies[2]*6])    
end


