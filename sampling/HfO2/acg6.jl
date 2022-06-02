cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf", "O"];
a = 3.59332800;
b = 3.59332800;
c = 5.22467300;
alpha = 90.0;
beta = 90;
gamma = 90.0;
atomtype = [1, 1, 2, 2, 2, 2];
fraccoords = [0.500000 0.500000 0.500000 
0.000000 0.000000 0.000000 
0.000000 0.500000 0.695308 
0.500000 0.000000 0.804692 
0.500000 0.000000 0.304692 
0.000000 0.500000 0.195308]';

lat = ACG.setlatticefraction(a, b, c, alpha, beta, gamma, fraccoords, atomtype);

rin = 0.5
rcut = 5.0
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Hf,:O], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Hf,:O], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Hf,:O], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors);

pranges = [a-0.7 a+0.7; b-0.7 b+0.7; c-1.0 c+1.0];        
pdims = [23, 23, 23];        
plattices, patombases = ACG.lengthconfigurations(lat, pranges, pdims, descriptors, 1500);

aranges = [75 105; 75 105; 75 105];        
adims = [23, 23, 23];        
alattices, aatombases = ACG.angleconfigurations(lat, aranges, adims, descriptors, 1500);

lranges = [a-0.7 a+0.7; b-0.7 b+0.7; c-1.0 c+1.0; 75 105; 75 105; 75 105];        
ldims = [6, 6, 6, 4, 4, 4];        
llattices, latombases = ACG.latticeconfigurations(lat, lranges, ldims, descriptors, 1500);

la = [plattices[:,:,1] alattices[:,:,1] llattices[:,:,1]]
ba = cat(patombases, aatombases, dims=3)
ba = cat(ba, latombases, dims=3)
dlattices, datombases = ACG.selectconfigurations(la[1:3,:], la[4:6,:], la[7:9,:], ba[:,:,:,1], atomtype, descriptors, 1e-6);

natom = length(atomtype);
ds = 0.075*ones(natom);        
N = 12000;    
blattices, batombases = ACG.atomdisplacementconfigurations(lat, ds, N, descriptors, 1e-5);
n = size(batombases,3)
blattices = blattices .* ones(9, n)

lattices = [dlattices[:,:,1] blattices[:,:,1]]
atombases = cat(datombases, batombases, dims=3)
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

ACG.writeEXTXYZ2("HfO2sg137n6.xyz", species, natoms, lattices, atomtypes, atombases);



