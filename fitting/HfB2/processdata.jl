cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

rin = 0.5
rcut = 5.0
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors)

# path to the database 
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Hf", "B"];

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

datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/"

folder = "EOS2";
folder = "Deformed3";
folder = "Distorted3"
folder = "MD2"
config = Preprocessing.readconfigpath(datapath * folder, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)

n = config.nconfigs;
m = config.natom[1];
la = config.lattice;
ba = reshape(config.x, (3, m, n));
atomtype = config.t[1:m];
lattices, atombases, ind, E = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba, atomtype, descriptors, 1e-8);

k = length(ind);
natom = config.natom[ind];
x = reshape(ba[:,:,ind], (3, m*k));
f = reshape(config.f, (3, m, n));
f = reshape(f[:,:,ind], (3, m*k));
e = config.e[ind]
s = reshape(config.stress[:,ind], (9, k))
t = reshape(config.t, (m, n))
t = reshape(t[:,ind], (1, m*k));

tconfig = initializeconfig(0);  
tconfig.natom = natom;
tconfig.x = x;
tconfig.t = t;
tconfig.f = f;
tconfig.e = e;
tconfig.stress = s;
tconfig.lattice = reshape(la[:,ind], (9, k))
tconfig.nconfigs = k;

filename = datapath * "training/" * folder * ".xyz";
Preprocessing.writeEXTXYZ(filename,tconfig,atomspecies);


datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/"

folder = "EOS2"
config = Preprocessing.readconfigpath(datapath * folder, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)
config.e = config.e * 27.211386;
config.f = config.f * 51.421;

filename = datapath * folder * "/" * folder * ".xyz";
Preprocessing.writeEXTXYZ(filename,config,atomspecies);


datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/HfB2_data_set/ACG/"
folder = "g6"
config = Preprocessing.readconfigpath(datapath * folder, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)
config.e = config.e * 27.211386;
config.f = config.f * 51.421;

datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/"
filename = datapath * "test" * "/" * folder * ".xyz";
Preprocessing.writeEXTXYZ(filename,config,atomspecies);


datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/HfB2_data_set/"
folder = "EOS2"
config = Preprocessing.readconfigpath(datapath * folder, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)
config.e = config.e * 27.211386;
config.f = config.f * 51.421;

n = config.nconfigs;
m = config.natom[1];
la = config.lattice;
ba = reshape(config.x, (3, m, n));
atomtype = config.t[1:m];
lattices, atombases, ind, E = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba, atomtype, descriptors, 1e-18);

k = length(ind);
natom = config.natom[ind];
x = reshape(ba[:,:,ind], (3, m*k));
f = reshape(config.f, (3, m, n));
f = reshape(f[:,:,ind], (3, m*k));
e = config.e[ind]
s = reshape(config.stress[:,ind], (9, k))
t = reshape(config.t, (m, n))
t = reshape(t[:,ind], (1, m*k));

tconfig = initializeconfig(0);  
tconfig.natom = natom;
tconfig.x = x;
tconfig.t = t;
tconfig.f = f;
tconfig.e = e;
tconfig.stress = s;
tconfig.lattice = reshape(la[:,ind], (9, k))
tconfig.nconfigs = k;

datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/"
filename = datapath * folder * "/" * folder * ".xyz";
Preprocessing.writeEXTXYZ(filename,tconfig,atomspecies);


datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/ACG/"
folder = "n81"
config = Preprocessing.readconfigpath(datapath * folder, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)

n = length(config.natom)
m = config.natom[1];
n1 = Int64(round(0.9*n));
for i = 1:2 
if (i==1)
    ind = 1:n1;
    datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/training/"
else
    ind = (n1+1):n;
    datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfB2/test/"
end
k = length(ind);
natom = config.natom[ind];
x = reshape(config.x, (3, m, n));
x = reshape(x[:,:,ind], (3, m*k));
f = reshape(config.f, (3, m, n));
f = reshape(f[:,:,ind], (3, m*k));
e = config.e[ind]
s = reshape(config.stress[:,ind], (9, k))
t = reshape(config.t, (m, n))
t = reshape(t[:,ind], (1, m*k));

tconfig = initializeconfig(0);  
tconfig.natom = reshape(natom,(1,k));
tconfig.x = x;
tconfig.t = t;
tconfig.f = f;
tconfig.e = reshape(e,(1,k));
tconfig.stress = s;
tconfig.lattice = reshape(config.lattice[:,ind], (9, k))
tconfig.nconfigs = k;

filename = datapath * "configs" * folder[2:end] * ".xyz";
Preprocessing.writeEXTXYZ(filename,tconfig,atomspecies);
end


# ACG.writeEXTXYZ2("Si/Si64.xyz", species, natoms, lattices, atomtypes, atombases);

# filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfO2/training/training.xyz";
# Preprocessing.writeEXTXYZ(filename,trainconfig,atomspecies);

# filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfO2/test/test.xyz";
# Preprocessing.writeEXTXYZ(filename,testconfig,atomspecies);

# filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfO2/validation/validation.xyz";
# Preprocessing.writeEXTXYZ(filename,validconfig,atomspecies);

# trainlat = config.lattice[:,ntrain]
# lat1 = unique(trainlat,dims=2);
# testlat = config.lattice[:,ntest]
# lat2 = unique(testlat,dims=2);





