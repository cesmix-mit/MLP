cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/HfO2/old_obsolete"
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Hf", "O"];

# weights for energies, forces, and stresses in the linear fit
weightinner=[100  1  1.00E-8]

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 1.0E-8]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = false 

config = Preprocessing.readconfigpath(datapath, dataformat, fileextension, atomspecies, transposelattice, nothing, 1)

lat = unique(config.lattice,dims=2);
o0 = [10.1973    9.55655   8.696    11.1311   10.0413   13.1104   11.7846     9.29098    8.34835  60.0
          0.0       0.0       0.0       0.0       0.0       0.0       0.0        0.0        0.0       0.0
          0.0       0.0       0.0       0.0       0.0       2.13638  -0.581793  -2.56984    1.60664   0.0
          0.0       0.0       8.696     2.78278   5.02064   7.85671  -2.2861     3.57928    5.56556   0.0
          8.34834  15.5452    7.90913  10.1973    8.696     8.74339  15.7527     7.32795   10.1973   60.0
          1.60664  -1.37401  -1.01934  -1.60664   3.09785  -1.01934  -0.99879    0.887048  -3.21327   0.0
          0.0       0.0       0.0       0.0       0.0       0.0       0.0        0.0        0.0       0.0
          0.0       0.0       0.0       0.0       0.0       0.0       0.0        0.0        0.0       0.0
          32.1328   29.5413   35.3967   32.1328   30.9784   35.3967   34.9712    36.1604    32.1328   60.0];

o1 = [11.1311   12.0131    12.0131     9.55655   9.55655    9.55655   9.01       9.55655   9.29353   9.08686   9.08686
0.0       0.0        0.0        0.0       0.0        0.0       0.0        0.0       0.0       0.0       0.0
0.0      -0.985032  -0.985032   0.0       0.0        0.0       0.0        0.0       0.0      -4.27275   4.27275
5.56556  -3.55448   -3.55448   -4.77827   0.0        0.0      -4.505      0.0       0.0       3.05958  -4.58937
10.1973   14.6559    14.6559     8.27621   9.55058    9.01     10.058      9.55655  13.044    10.0919   15.1379
-3.21328  -0.676355  -0.676355   0.0      -0.337661   0.0      -0.558773   0.0       2.51032  -5.29209  -2.03868
0.0       0.0        0.0        0.0       0.0        0.0       0.0        0.0       0.0       0.0       0.0
0.0       0.0        0.0        0.0       0.0        0.0       0.0        0.0       0.0       0.0       0.0
32.1328   31.0739    31.0739    31.2115   30.0521    38.2262   36.3205    31.535    35.1445   35.3967   35.3967]          

o2 = [8.696     8.696     8.696     8.696    10.0413   10.0413   10.0413    8.696     9.29353  11.1311   10.1973
0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
-5.02064   5.02064   0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
0.0       0.0       0.0       8.696    -5.02064  -5.02064   5.02064  -8.696     0.0       0.0       0.0
9.29354   9.29354  10.5455   13.1819    8.696     8.696     8.696    13.1819    8.69601  10.8738   11.1311
0.0       0.0       5.29208   2.23406   3.09785   0.0       0.0       2.23406   0.0      -5.20055   0.0
0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
35.1445   35.1445   35.3967   35.3967   30.9784   30.9784   30.9784   35.3967   30.1238   36.1604   32.1328]

outlat = [o0 o1 o2]


ind = zeros(Int32,size(lat,2))
for i = 1:size(lat,2)
    for j = 1:size(outlat,2)        
        tm = outlat[:,j] - lat[:,i]
        a = sum(tm.^2);
        if a<1e-4
            ind[i] = 1
            break;
        end
    end    
end
ind = findall(ind .== 0)
lat = lat[:,ind]

nlt = size(lat,2)
ind = Array{Any}(nothing, nlt)
n = zeros(Int32, nlt)
for i = 1:nlt
    tm = config.lattice .- lat[:,i]
    a = sum(tm.^2,dims=1);
    ind[i] = findall(a[:] .< 1e-12)
    n[i] = length(ind[i])
end

training = 0.04;
test = 0.008; 
valid = 0.008; 

ntrain = []
ntest = []
nvalid = []
for i = 1:nlt
    n1 = Int32(round(n[i]*training)) 
    n2 = Int32(round(n[i]*test)) 
    n3 = Int32(round(n[i]*valid)) 
    a1 = randperm(n[i])
    a2 = randperm(n[i])
    a3 = randperm(n[i])
    b1 = a1[1:n1]
    b2 = a2[1:n2]
    b3 = a3[1:n3]
    ntrain = [ntrain; ind[i][b1]]
    ntest = [ntest; ind[i][b2]]
    nvalid = [nvalid; ind[i][b3]]
end

nat = [0; cumsum(config.natom[:])];

trainconfig = initializeconfig(0);  
testconfig = initializeconfig(0);  
validconfig = initializeconfig(0);  

trainconfig.nconfigs = length(ntrain)
trainconfig.natom = reshape(config.natom[ntrain], (1, trainconfig.nconfigs))
trainconfig.lattice = config.lattice[:,ntrain]
trainconfig.e = reshape(config.e[ntrain], (1, trainconfig.nconfigs))
trainconfig.stress = config.stress[:,ntrain]
for j = 1:length(ntrain)
    i = ntrain[j]
    t = config.t[:,(nat[i]+1):nat[i+1]]
    x = config.x[:,(nat[i]+1):nat[i+1]]
    f = config.f[:,(nat[i]+1):nat[i+1]]
    if j == 1
        trainconfig.t = t
        trainconfig.x = x
        trainconfig.f = f        
    else
        trainconfig.t = [trainconfig.t t]
        trainconfig.x = [trainconfig.x x]
        trainconfig.f = [trainconfig.f f]
    end
    d = trainconfig.lattice[:,j] - config.lattice[:,i];
    if sum(d.^2) > 1e-10
        error("something wrong")
    end
end

testconfig.nconfigs = length(ntest)
testconfig.natom = reshape(config.natom[ntest], (1, testconfig.nconfigs))
testconfig.lattice = config.lattice[:,ntest]
testconfig.e = reshape(config.e[ntest], (1, testconfig.nconfigs))
testconfig.stress = config.stress[:,ntest]
for j = 1:length(ntest)
    i = ntest[j]
    t = config.t[:,(nat[i]+1):nat[i+1]]
    x = config.x[:,(nat[i]+1):nat[i+1]]
    f = config.f[:,(nat[i]+1):nat[i+1]]
    if j == 1
        testconfig.t = t
        testconfig.x = x
        testconfig.f = f        
    else
        testconfig.t = [testconfig.t t]
        testconfig.x = [testconfig.x x]
        testconfig.f = [testconfig.f f]
    end
    d = testconfig.lattice[:,j] - config.lattice[:,i];
    if sum(d.^2) > 1e-10
        error("something wrong")
    end
end

validconfig.nconfigs = length(nvalid)
validconfig.natom = reshape(config.natom[nvalid], (1, validconfig.nconfigs))
validconfig.lattice = config.lattice[:,nvalid]
validconfig.e = reshape(config.e[nvalid], (1, validconfig.nconfigs))
validconfig.stress = config.stress[:,nvalid]
for j = 1:length(nvalid)
    i = nvalid[j]
    t = config.t[:,(nat[i]+1):nat[i+1]]
    x = config.x[:,(nat[i]+1):nat[i+1]]
    f = config.f[:,(nat[i]+1):nat[i+1]]
    if j == 1
        validconfig.t = t
        validconfig.x = x
        validconfig.f = f        
    else
        validconfig.t = [validconfig.t t]
        validconfig.x = [validconfig.x x]
        validconfig.f = [validconfig.f f]
    end
    d = validconfig.lattice[:,j] - config.lattice[:,i];
    if sum(d.^2) > 1e-10
        error("something wrong")
    end
end

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfO2/training/training.xyz";
Preprocessing.writeEXTXYZ(filename,trainconfig,atomspecies);

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfO2/test/test.xyz";
Preprocessing.writeEXTXYZ(filename,testconfig,atomspecies);

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/HfO2/validation/validation.xyz";
Preprocessing.writeEXTXYZ(filename,validconfig,atomspecies);


# trainlat = config.lattice[:,ntrain]
# lat1 = unique(trainlat,dims=2);
# testlat = config.lattice[:,ntest]
# lat2 = unique(testlat,dims=2);





