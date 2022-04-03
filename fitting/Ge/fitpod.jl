cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/Ge/"
dataformat = "jsonpymatgen"
fileextension = "json"
atomspecies = ["Ge"];

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

# training data 
traindata[1] = adddata(datapath * "training", dataformat, fileextension, 
            percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
            rotationmatrix, transposelattice)

# testing data 
testdata[1] = adddata(datapath * "test", dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[1,:], translationvector, 
                rotationmatrix, transposelattice)

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

eta = zeros(6,2)
for j = 1:6
    tm = readdlm("results/optpodeta" * string(j-1) *  ".txt")
    eta[j,:] = tm[:]
end
#show(IOContext(stdout, :compact=>false), "text/plain", eta)

# Bessel scaling parameters
gamma = [0.0, 2, 4]
for j = 0:5    

    rcut = eta[j+1,1]
    rin = eta[j+1,2]
    display([j rcut rin])
    
    # POD Descriptors    
    if j == 0
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j == 1                
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==2
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==3 
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==4
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==5        
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)   
        descriptors[2] = POD(nbody=2, pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[3] = POD(nbody=3, pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma)
    end
    
    coeff = linearfit(traindata, descriptors, Doptions)
    energyerrors, forceerrors = Potential.poderroranalysis(traindata, descriptors, Doptions, coeff)
    energytesterrors, forcetesterrors = Potential.poderroranalysis(testdata, descriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(["test"], e1, "Energy Errors")
    printerrors(["test"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitpodcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitpodtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/fitpodtesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
end


