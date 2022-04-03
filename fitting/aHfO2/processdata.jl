cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Potential Comparison/data/a-HfO2"
#datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/HfO2/old_obsolete"
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

nlt = size(lat,2)
ind = Array{Any}(nothing, nlt)
n = zeros(Int32, nlt)
for i = 1:nlt
    tm = config.lattice .- lat[:,i]
    a = sum(tm.^2,dims=1);
    ind[i] = findall(a[:] .< 1e-12)
    n[i] = length(ind[i])
end

training = 0.25;
test = 0.05; 
valid = 0.05; 

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

function forcemean(f, natom)
    nat = [0; cumsum(natom[:])];
    n = length(natom)
    fm = 0;
    for i = 1:n        
        fm = fm .+ f[:,(nat[i]+1):nat[i+1]]        
    end
    fm = fm/n;
    return fm
    # fn = 0.0*f;
    # for i = 1:n
    #     fn[:,(nat[i]+1):nat[i+1]] = f[:,(nat[i]+1):nat[i+1]] - fm         
    # end
    # return fm, fn 
end

function subtractmean(f, fm, natom)
    nat = [0; cumsum(natom[:])];
    n = length(natom)
    fn = 0.0*f;
    for i = 1:n
        fn[:,(nat[i]+1):nat[i+1]] = f[:,(nat[i]+1):nat[i+1]] - fm         
    end
    return fn
end

fm = forcemean(config.f, config.natom)
em = sum(config.e)/length(config.natom) 

trainconfig.e = trainconfig.e .- em; 
testconfig.e = testconfig.e .- em; 
validconfig.e = validconfig.e .- em; 

trainconfig.f = subtractmean(trainconfig.f, fm, trainconfig.natom);
testconfig.f = subtractmean(testconfig.f, fm, testconfig.natom);
validconfig.f = subtractmean(validconfig.f, fm, validconfig.natom);

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/aHfO2/training/training.xyz";
Preprocessing.writeEXTXYZ(filename,trainconfig,atomspecies);

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/aHfO2/test/test.xyz";
Preprocessing.writeEXTXYZ(filename,testconfig,atomspecies);

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MLP/data/aHfO2/validation/validation.xyz";
Preprocessing.writeEXTXYZ(filename,validconfig,atomspecies);

# # trainlat = config.lattice[:,ntrain]
# # lat1 = unique(trainlat,dims=2);
# # testlat = config.lattice[:,ntest]
# # lat2 = unique(testlat,dims=2);





