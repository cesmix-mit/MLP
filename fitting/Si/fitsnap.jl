cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

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

twojmaxlist = [2, 4, 6, 8, 10, 12]
for j = 0:5        
    display(j)
    twojmx = twojmaxlist[j+1]
    
    if j == 0
        rcut = 4.5                
    elseif j == 1
        rcut = 4.5        
    elseif j==2
        rcut = 5.276384856307164  
    elseif j==3 
        rcut = 5.227869204194917   
    elseif j==4
        rcut = 5.23859140860721
    elseif j==5           
        rcut = 5.234331023785782
    end
    descriptors[1] = SNAPparams(species = [:Si], nbody = 1, twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363)
    descriptors[2] = SNAPparams(species = [:Si], nbody = 4, twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363)
    
    coeff = linearfit(traindata, descriptors, Doptions)
    energyerrors, forceerrors = Potential.snaperroranalysis(traindata, descriptors, Doptions, coeff)
    energytesterrors, forcetesterrors = Potential.snaperroranalysis(testdata, descriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(["test"], e1, "Energy Errors")
    printerrors(["test"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitsnapcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitsnaptrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/fitsnaptesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
end

# eta = zeros(6,1)
# for j = 1:6
#     tm = readdlm("results/optsnapeta" * string(j-1) *  ".txt")
#     eta[j] = tm[1]
# end
# show(IOContext(stdout, :compact=>false), "text/plain", eta)

