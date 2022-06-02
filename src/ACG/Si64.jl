cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Si"];
m = 2;
n = 2;
p = 2;
a = 5.43*m;
b = 5.43*n;
c = 5.43*p;
alpha = 90;
beta = 90;
gamma = 90;
fraccoords = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5; 0.25 0.25 0.25; 0.25 0.75 0.75; 0.75 0.25 0.75; 0.75 0.75 0.25]
fraccoords = ACG.replicatefraccoords(fraccoords, m, n, p);
atomtype = ones(Int64, 64) 

lat = ACG.setlatticefraction(a, b, c, alpha, beta, gamma, fraccoords, atomtype);

rin = 1.0
rcut = 5.2
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Si], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Si], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Si], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors)

pranges = [6.867 15.867; 6.867 15.867; 6.867 15.867];        
pdims = [21, 21, 21];        
plattices, patombases = ACG.lengthconfigurations(lat, pranges, pdims, descriptors);

aranges = [60 115; 60 115; 60 90];        
adims = [21, 21, 21];        
alattices, aatombases = ACG.angleconfigurations(lat, aranges, adims, descriptors);

lranges = [6.867 15.867; 6.867 15.867; 6.867 15.867; 60 115; 60 115; 60 90];        
ldims = [6, 6, 6, 4, 4, 4];        
llattices, latombases = ACG.latticeconfigurations(lat, lranges, ldims, descriptors);

ds = 0.05*ones(64);        
N = 10000;    
blattices, batombases = ACG.atomdisplacementconfigurations(lat, ds, N, descriptors);
n = size(batombases,3)
blattices = blattices .* ones(9, n)

la = [plattices[:,:,1] alattices[:,:,1] llattices[:,:,1] blattices[:,:,1]]
ba = cat(patombases, aatombases, dims=3)
ba = cat(ba, latombases, dims=3)
ba = cat(ba, batombases, dims=3)

lattices, atombases = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba[:,:,:,1], atomtype, descriptors);

n = size(lattices, 2);
natoms = 64*ones(Int64, n);
lattices = reshape(lattices, (9, n))
atombases = reshape(atombases, (3, 64*n));
atomtypes = ones(Int64, 64*n)

ACG.writeEXTXYZ2("Si/Si64.xyz", species, natoms, lattices, atomtypes, atombases);

