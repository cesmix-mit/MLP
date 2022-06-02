cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf", "O"];
a = 5.18420400;
b = 5.29677200;
c = 10.16892500;
alpha = 90.0;
beta = 90.0;
gamma = 90.0;
atomtype = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
fraccoords = [0.958737 0.342789 0.137946 
0.541263 0.842789 0.137946 
0.458737 0.342789 0.362054 
0.041263 0.842789 0.362054 
0.958737 0.157211 0.637946 
0.541263 0.657211 0.637946 
0.458737 0.157211 0.862054 
0.041263 0.657211 0.862054 
0.832426 0.667500 0.033559 
0.667574 0.167500 0.033559 
0.248913 0.087304 0.224333 
0.251087 0.587304 0.224333 
0.748913 0.087304 0.275667 
0.751087 0.587304 0.275667 
0.167574 0.167500 0.466441 
0.332426 0.667500 0.466441 
0.667574 0.332500 0.533559 
0.832426 0.832500 0.533559 
0.248913 0.412696 0.724333 
0.251087 0.912696 0.724333 
0.748913 0.412696 0.775667 
0.751087 0.912696 0.775667 
0.167574 0.332500 0.966441 
0.332426 0.832500 0.966441]';

ranges = [a-1.0 a+1.0; b-1.0 b+1.0];        
dims = [25, 25];        
s = Optimization.tensornodes(ranges[:,1], ranges[:,2], dims)

natom = length(atomtype);
n = size(s, 1);
lattices = zeros(9,n)
atombases = zeros(3, natom, n)
for i = 1:n
    lat = ACG.setlatticefraction(s[i,1], s[i,2], c, alpha, beta, gamma, fraccoords, atomtype);
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

ACG.writeEXTXYZ2("HfO2sg61n24test.xyz", species, natoms, lattices, atomtypes, atombases);

