cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Si"];
a = 3.867;
b = 3.867;
c = 3.867;
alpha = 60;
beta = 60;
gamma = 60;
atomtype = [1, 1]
fraccoords = [0.0  0.0  0.0; 
             0.25 0.25 0.25]

lat = ACG.setlatticefraction(a, b, c, alpha, beta, gamma, fraccoords, atomtype);

rin = 1.0
rcut = 5.2
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Si], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Si], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Si], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors)

pranges = [2.367 5.367; 2.367 5.367; 2.367 5.367];        
pdims = [21, 21, 21];        
plattices, patombases = ACG.lengthconfigurations(lat, pranges, pdims, descriptors);

aranges = [45 80; 45 80; 45 80];        
adims = [21, 21, 21];        
alattices, aatombases = ACG.angleconfigurations(lat, aranges, adims, descriptors);

lranges = [2.367 5.367; 2.367 5.367; 2.367 5.367; 45 80; 45 80; 45 80];        
ldims = [6, 6, 6, 4, 4, 4];        
llattices, latombases = ACG.latticeconfigurations(lat, lranges, ldims, descriptors);

ds = 0.1*ones(2);        
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
natoms = 2*ones(Int64, n);
lattices = reshape(lattices, (9, n))
atombases = reshape(atombases, (3, 2*n));
atomtypes = ones(Int64, 2*n)

ACG.writeEXTXYZ2("Si/Si2.xyz", species, natoms, lattices, atomtypes, atombases);


# branges = [0.0, 0.125];        
# bdims = [1, 16];        
# blattices, batombases = ACG.atombasisconfigurations(lat, branges, bdims, descriptors);
# n = size(batombases,3)
# blattices = blattices .* ones(9, n)

# j = 108;
# a1 = lattices[1:3,j]
# a2 = lattices[4:6,j]
# a3 = lattices[7:9,j]
# pos = [a1 a2 a3]*fraccoords'
# atombases[:,:,j]-pos
# Preprocessing.mkfolder("Si")
# writedlm("Si/lattices.txt", lattices[:,:,1]')





