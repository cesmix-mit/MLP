cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

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

eta = zeros(6,1)
for j = 1:6
    tm = readdlm("results/optsnapeta" * string(j-1) *  ".txt")
    eta[j] = tm[1]
end
# show(IOContext(stdout, :compact=>false), "text/plain", eta)

twojmaxlist = [2, 4, 6, 8, 10, 12]
for j = 0:5        
    twojmx = twojmaxlist[j+1]    
    rcut = eta[j+1]
    display([j twojmx rcut])

    descriptors[1] = SNAPparams(species = [:Ta], nbody = 1, twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363)
    descriptors[2] = SNAPparams(species = [:Ta], nbody = 4, twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363)
    
    coeff = linearfit(traindata, descriptors, Doptions)
    energyerrors, forceerrors, stresserrors, eerr, ferr, ferm, energies = Potential.snaperroranalysis(traindata, descriptors, Doptions, coeff)

    show(stdout, "text/plain", energies[11])

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(folders, e1, "Energy Errors")
    printerrors(folders, e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitsnapcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitsnaptrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
end

