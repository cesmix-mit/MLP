cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf", "O"];
a = 4.93704200;
b = 5.13372400;
c = 5.21084800;
alpha = 90.0;
beta = 90.0;
gamma = 90.0;
atomtype = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2];
fraccoords = [0.728074 0.277502 0.282356 
0.271926 0.777502 0.217644 
0.728074 0.722498 0.717644 
0.271926 0.222498 0.782356 
0.425916 0.500000 0.500000 
0.574084 0.000000 0.000000 
0.510630 0.500000 0.000000 
0.489370 0.000000 0.500000 
0.896826 0.313870 0.662984 
0.103174 0.813870 0.837016 
0.896826 0.686130 0.337016 
0.103174 0.186130 0.162984]';

ranges = [a-1.0 a+1.0; c-1.0 c+1.0];        
dims = [25, 25];        
s = Optimization.tensornodes(ranges[:,1], ranges[:,2], dims)

natom = length(atomtype);
n = size(s, 1);
lattices = zeros(9,n)
atombases = zeros(3, natom, n)
for i = 1:n
    lat = ACG.setlatticefraction(s[i,1], b, s[i,2], alpha, beta, gamma, fraccoords, atomtype);
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

ACG.writeEXTXYZ2("HfO2sg18n12test.xyz", species, natoms, lattices, atomtypes, atombases);

