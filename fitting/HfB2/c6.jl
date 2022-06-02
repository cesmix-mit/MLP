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

ranges = [a-0.5 a+0.5; c-0.5 c+0.5];        
dims = [24, 24];        
s = Optimization.tensornodes(ranges[:,1], ranges[:,2], dims)

natom = length(atomtype);
n = size(s, 1);
lattices = zeros(9,n)
atombases = zeros(3, natom, n)
for i = 1:n
    lat = ACG.setlatticefraction(s[i,1], s[i,1]*sqrt(3), s[i,2], alpha, beta, gamma, fraccoords, atomtype);
    lattices[:,i] = [lat.a1; lat.a2; lat.a3]
    atombases[:,:,i] = lat.atombasis
end

natoms = natom*ones(Int64, n);
lattices = reshape(lattices, (9, n))
atombases = reshape(atombases, (3, natom*n));
atomtypes = atomtype[:] .* ones(Int64, (1, n)); 
atomtypes = reshape(atomtypes, (1, natom*n))

ind = findall(abs.(lattices[:]) .< 1e-12);
lattices[ind] .= 0.0;
ind = findall(abs.(atombases[:]) .< 1e-12);
atombases[ind] .= 0.0;

ACG.writeEXTXYZ2("HfB2g6.xyz", species, natoms, lattices, atomtypes, atombases);

