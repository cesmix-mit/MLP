cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf"];
# path to the database 
dataformat = "extxyz"
datamode = 1
fileextension = "xyz"
atomspecies = ["Hf"];

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = false 

datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/HfB2_data_set/HfO2/Slabs";

datafile = "Hfmp103n128.xyz";
config = readconfigfile(datapath, datafile, dataformat, datamode, atomspecies, transposelattice, nothing, 1);   
atomtype = ones(Int64, (1, 128)); 
la = config.lattice;
n = length(config.natom)
ba = reshape(config.x, (3, 128, n));

rin = 1.0
rcut = 5.0
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Hf], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Hf], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Hf], pdegree=[3,6,3], nbasis = [4, 3], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors);

lattices, atombases = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba, atomtype, descriptors, 2e-4);

natom = 128;
n = size(lattices, 2);
natoms = natom*ones(Int64, n);
lattices = reshape(lattices, (9, n))
atombases = reshape(atombases, (3, natom*n));
atomtypes = atomtype[:] .* ones(Int64, (1, n)); 
atomtypes = reshape(atomtypes, (1, natom*n))

ind = findall(abs.(lattices[:]) .< 1e-12);
lattices[ind] .= 0.0;
ind = findall(abs.(atombases[:]) .< 1e-12);
atombases[ind] .= 0.0;

ACG.writeEXTXYZ2("Hfmp103n128_selected.xyz", species, natoms, lattices, atomtypes, atombases);

# end 