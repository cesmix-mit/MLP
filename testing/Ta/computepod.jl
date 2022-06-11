cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/Ta/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["Ta"];

# folders in the datapath 
# folders = ["Displaced_A15", "Displaced_BCC", "Displaced_FCC", "Elastic_BCC",
#             "Elastic_FCC", "GSF_110", "GSF_112", "Liquid", "Surface",
#             "Volume_A15", "Volume_BCC", "Volume_FCC"];

folders = ["Displaced_A15"];

# weights for energies, forces, and stresses in the linear fit
weightinner=[100            1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09
            100             1               1.00E-09]

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 1.0E-9]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = true 

# training data 
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
                rotationmatrix, transposelattice)
end

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# optimization parameters
optim = setoptim(lossfunc, method)

# Bessel scaling parameters
gamma = [0.0, 2, 4]
for j = 2:2
    
    # POD Descriptors
    if j == 0
        rcut = 4.653843683095952
        rin = 1.0405983377652281
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j == 1
        rcut = 4.4902344770267195
        rin = 1.4229278161004044
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==2
        rcut = 5.0
        rin = 1.145304772214426
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==3 
        rcut = 4.6829290489590045
        rin = 0.9
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==4
        rcut = 4.748396818189743
        rin = 0.9
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==5        
        rcut = 4.843466804257663
        rin = 0.9364496626725823
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)   
        descriptors[2] = POD(nbody=2, pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma)
    end
    
    Preprocessing.mkfolder("results")
    for i = 1:length(folders)        
        config, ~, ~ = Preprocessing.readconfigdata(traindata[i], Doptions.pbc)    
        natom = [0; cumsum(config.natom[:])];
        for n = 1:length(config.natom)
            x, t, a, b, c, pbc, e, f = Potential.getconfig(config, natom, n)
            dglob, dderiv, datom = Potential.podefatom(x, t, a, b, c, pbc, descriptors)          
            writedlm("results/computepoddglob" * string(j) * string(i) * string(n) *  ".txt", dglob)
            writedlm("results/computepoddderiv" * string(j) * string(i) * string(n) *  ".txt", dderiv)
            writedlm("results/computepoddatom" * string(j) * string(i) * string(n) *  ".txt", datom) 
        end
    end
end

