cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles
using ACG, Potential

species = ["Hf"];
a = 3.19848594;
b = 3.19848594;
c = 5.07518500;
alpha = 90.0;
beta = 90;
gamma = 120.0;
atomtype = [1, 1];
fraccoords = [0.666667 0.333333 0.250000
0.333333 0.666667 0.750000]';

lat = ACG.setlatticefraction(a, b, c, alpha, beta, gamma, fraccoords, atomtype);

rin = 0.5
rcut = 5.0
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Hf], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Hf], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Hf], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors);

pranges = [a-0.7 a+0.7; b-0.7 b+0.7; c-1.0 c+1.0];        
pdims = [23, 23, 23].-1;        

aranges = [75.0 105.0; 75 105; 100 140]*pi/180.0;        
adims = [23, 23, 23].-1;        

lranges = [a-0.7 a+0.7; b-0.7 b+0.7; c-1.0 c+1.0; aranges];        
ldims = [6, 6, 6, 4, 4, 4].-1;        

natom = length(atomtype);
ds = 0.075*ones(natom);        
N = 15000;    

Preprocessing.mkfolder("results")

A1, A2, A3, Atombases = ACG.lengthgroup(lat, pranges, pdims);
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
savebin("Hfg194n6e_lengthgroup.bin", edesc);
savebin("Hfg194n6g_lengthgroup.bin", gdesc);

A1, A2, A3, Atombases = ACG.anglegroup(lat, aranges, adims);        
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
savebin("Hfg194n6e_anglegroup.bin", edesc);
savebin("Hfg194n6g_anglegroup.bin", gdesc);

A1, A2, A3, Atombases = ACG.latticegroup(lat, lranges, ldims);
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
savebin("Hfg194n6e_latticegroup.bin", edesc);
savebin("Hfg194n6g_latticegroup.bin", gdesc);

A1, A2, A3, Atombases = ACG.atomdisplacementgroup(lat, ds, N);
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
savebin("Hfg194n6e_atomgroup.bin", edesc);
savebin("Hfg194n6g_atomgroup.bin", gdesc);

