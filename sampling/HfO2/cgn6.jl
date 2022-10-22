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
pdims = [23, 23, 23].-1;        

aranges = [75.0 105.0; 75 105; 75 105]*pi/180.0;        
adims = [23, 23, 23].-1;        

lranges = [a-0.7 a+0.7; b-0.7 b+0.7; c-1.0 c+1.0; aranges];        
ldims = [6, 6, 6, 4, 4, 4].-1;        

natom = length(atomtype);
ds = 0.075*ones(natom);        
N = 15000;    

Preprocessing.mkfolder("results")

A1, A2, A3, Atombases = ACG.lengthgroup(lat, pranges, pdims);
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
#writedlm("results/HfO2g137n6e_lengthgroup.txt", edesc);
writedlm("results/HfO2g137n6g_lengthgroup.txt", gdesc);

A1, A2, A3, Atombases = ACG.anglegroup(lat, aranges, adims);        
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
#writedlm("results/HfO2g137n6e_anglegroup.txt", edesc);
writedlm("results/HfO2g137n6g_anglegroup.txt", gdesc);

A1, A2, A3, Atombases = ACG.latticegroup(lat, lranges, ldims);
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
#writedlm("results/HfO2g137n6e_latticegroup.txt", edesc);
writedlm("results/HfO2g137n6g_latticegroup.txt", gdesc);

A1, A2, A3, Atombases = ACG.atomdisplacementgroup(lat, ds, N);
edesc, gdesc  = ACG.latticepoddescriptors(A1, A2, A3, Atombases, atomtype, descriptors)
#writedlm("results/HfO2g137n6e_atomgroup.txt", edesc);
writedlm("results/HfO2g137n6g_atomgroup.txt", gdesc);

