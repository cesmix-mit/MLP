cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/TiO2"
dataformat = "ann"
fileextension = "xsf"
atomspecies = ["Ti", "O"];

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

config = Preprocessing.readconfigpath(datapath, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)
config.stress = zeros(9, config.nconfigs);

Preprocessing.writeEXTXYZ("TiO2.xyz",config,atomspecies);

percentage=50;
randomize=0;
config.nconfigs = length(config.natom)
ntrain = Int64(round(percentage*config.nconfigs/100.0)) 
indices = Array(1:ntrain)     
if randomize==1
    ind = randperm(config.nconfigs)
    indices = ind[1:ntrain]        
elseif randomize==0
    indices = Int32.(round.(LinRange(1, config.nconfigs, ntrain)))            
end          
trainconfig = Preprocessing.extractconfig(config, indices);
Preprocessing.writeEXTXYZ("TiO2trainingset.xyz",trainconfig,atomspecies);


ind = randperm(config.nconfigs);
indices = ind[1:ntrain]        
testconfig = Preprocessing.extractconfig(config, indices);
Preprocessing.writeEXTXYZ("TiO2testset.xyz",testconfig,atomspecies);



# savebin("natom.bin",config.natom[:])
# savebin("x.bin",config.x[:])
# savebin("t.bin",config.t[:])
# savebin("e.bin",config.e[:])
# savebin("lattice.bin",config.lattice[:])

