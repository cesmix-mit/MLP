cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

using ACG, Potential

pbc = [1, 1, 1]
lat = [2.96588026 0.0 0.0 0.0 4.6471346 0.0 0.0 -4.022e-5 4.64713459]
atomspecies = ["Ti", "O"]
atomtype = [1, 1, 2, 2, 2, 2]
atombasis = [1.48294008  2.32354719  2.3235673; 
             0.0  0.0  0.0; 
             1.48294008  0.90569453  0.90570237;
             1.48294008  3.74139995  3.74143233;
             0.0  1.41783138  3.22927527;
             0.0  3.229263  1.41785932]

a1 = lat[1:3]
a2 = lat[4:6]
a3 = lat[7:9]

lat = ACG.setlattice("custom", a1, a2, a3, atombasis, atomtype);

rin = 1.0
rcut = 5.2
gamma = [0.0, 2, 4]
descriptors[1] = POD(nbody=1, species = [:Ti,:O], pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
descriptors[2] = POD(nbody=2, species = [:Ti,:O], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors[3] = POD(nbody=3, species = [:Ti,:O], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
descriptors = Preprocessing.deletenothing(descriptors)

ranges = [2.0 4.0; 3.75 5.5; 3.75 5.5];        
dims = [11, 11, 11].-1;        
A1, A2, A3, AtomBase = ACG.lengthgroup(lat, ranges, dims);
E, F, G, edesc, fdesc, gdesc  = ACG.latticesimilaritymatrices(A1, A2, A3, AtomBase, lat.atomtype, descriptors)
indE = ACG.minmaxsampling(E, 100);
indF = ACG.minmaxsampling(F, 100);
indG = ACG.minmaxsampling(G, 100);
ind = ACG.minmaxsampling(G.*F, 100);


ranges = [75 105; 75 105; 75 105]*pi/180.0;        
dims = [25, 25, 25].-1;        
A1, A2, A3, AtomBase = ACG.anglegroup(lat, ranges, dims);

ranges = [1.7 3.9; 3.0 5.9; 3.0 5.9; 80 100; 80 100; 80 100]*pi/180.0;        
dims = [20, 20, 20, 4, 4, 4].-1;        
A1, A2, A3, AtomBase = ACG.latticegroup(lat, ranges, dims);

nbasis = size(lat.fraccoords,2)
ds = 0.1*ones(Float64,nbasis);     

dims = [3, 3, 1, 1, 1, 1];        
A1, A2, A3, AtomBase = ACG.atombasisgroup(lat, ds, dims);

dims = [1, 1, 3, 3, 1, 1];        
A1, A2, A3, AtomBase = ACG.atombasisgroup(lat, ds, dims);

dims = [1, 1, 1, 1, 3, 3];        
A1, A2, A3, AtomBase = ACG.atombasisgroup(lat, ds, dims);


# Lattice="2.96588026 0.0 0.0 0.0 4.6471346 0.0 0.0 -4.022e-5 4.64713459" Properties=species:S:1:pos:R:3:forces:R:3 energy=-4990.44928914 stress="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" pbc="T T T"
# Ti  1.48294008  2.32354719  2.3235673  0.0  0.0  0.0
# Ti  0.0  0.0  0.0  0.0  0.0  0.0
# O  1.48294008  0.90569453  0.90570237  -0.0  -0.00170978  -0.00170978
# O  1.48294008  3.74139995  3.74143233  0.0  0.00170978  0.00170978
# O  0.0  1.41783138  3.22927527  0.0  0.00197872  -0.00197872
# O  0.0  3.229263  1.41785932  0.0  -0.00197872  0.00197872



