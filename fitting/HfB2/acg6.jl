cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf", "B"];
a = 3.149;
b = 3.149*sqrt(3);
c = 3.480;
alpha = 90;
beta = 90;
gamma = 90;
atomtype = [1, 1, 2, 2, 2, 2];
fraccoords = [0 0 0; 1/2 1/2 0; 1/2 1/6 1/2; 0 2/6 1/2; 0 4/6 1/2; 1/2 5/6 1/2]';

lat = ACG.setlatticefraction(a, b, c, alpha, beta, gamma, fraccoords, atomtype);

rin = 1.0
rcut = 5.0
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Hf,:B], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Hf,:B], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Hf,:B], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors);

pranges = [a-0.5 a+0.5; b-0.8 b+0.8; c-0.5 c+0.5];        
pdims = [21, 21, 21];        
plattices, patombases = ACG.lengthconfigurations(lat, pranges, pdims, descriptors, 1e-6);

aranges = [70 90; 80 100; 80 100];        
adims = [21, 21, 21];        
alattices, aatombases = ACG.angleconfigurations(lat, aranges, adims, descriptors, 1e-6);

lranges = [a-0.5 a+0.5; b-0.8 b+0.8; c-0.5 c+0.5; 70 90; 80 100; 80 100];        
ldims = [6, 6, 6, 4, 4, 4];        
llattices, latombases = ACG.latticeconfigurations(lat, lranges, ldims, descriptors, 1e-6);

natom = length(atomtype);
ds = 0.05*ones(natom);        
N = 10000;    
blattices, batombases = ACG.atomdisplacementconfigurations(lat, ds, N, descriptors, 1e-6);
n = size(batombases,3)
blattices = blattices .* ones(9, n)

la = [plattices[:,:,1] alattices[:,:,1] llattices[:,:,1] blattices[:,:,1]]
ba = cat(patombases, aatombases, dims=3)
ba = cat(ba, latombases, dims=3)
ba = cat(ba, batombases, dims=3)

lattices, atombases = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba[:,:,:,1], atomtype, descriptors, 5e-5);

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

ACG.writeEXTXYZ2("HfB2n6.xyz", species, natoms, lattices, atomtypes, atombases);


