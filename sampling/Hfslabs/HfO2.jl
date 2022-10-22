cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf", "O"];
# path to the database 
dataformat = "extxyz"
datamode = 1
fileextension = "xyz"
atomspecies = ["Hf", "O"];

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = false 

datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/HfB2_data_set/HfO2/Slabs";

datafile = "HfO2mp1018721n162.xyz";
config = readconfigfile(datapath, datafile, dataformat, datamode, atomspecies, transposelattice, nothing, 1);   
#atomtype = ones(Int64, (1, 250)); 
la = config.lattice;
n = length(config.natom)
t = reshape(config.t, (config.natom[1], n));
atomtype = t[:,1]
ba = reshape(config.x, (3, config.natom[1], n));

rin = 1.0
rcut = 5.0
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Hf, :O], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Hf, :O], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Hf, :O], pdegree=[3,6,3], nbasis = [4, 3], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors);

lattices, atombases = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba, atomtype, descriptors, 2e-4);

natom = config.natom[1];
n = size(lattices, 2);
atombases = reshape(atombases, (3, natom*n))
natoms = natom*ones(Int64, n);
lattices = reshape(lattices[:,1:n], (9, n))
atombases = reshape(atombases[:,1:natom*n], (3, natom*n));
atomtypes = atomtype[:] .* ones(Int64, (1, n)); 
atomtypes = reshape(atomtypes, (1, natom*n))

ind = findall(abs.(lattices[:]) .< 1e-12);
lattices[ind] .= 0.0;
ind = findall(abs.(atombases[:]) .< 1e-12);
atombases[ind] .= 0.0;

ACG.writeEXTXYZ2("HfO2mp1018721n162_selected.xyz", species, natoms, lattices, atomtypes, atombases);

