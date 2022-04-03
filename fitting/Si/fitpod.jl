cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/Si/"
dataformat = "jsonpymatgen"
fileextension = "json"
atomspecies = ["Si"];

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

# Bessel scaling parameters
gamma = [0.0, 2, 4]
for j = 0:5    
    # POD Descriptors
    if j == 0
        rcut = 4.279154036156406
        rin = 0.9
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
    elseif j == 1
        rcut = 4.287842639403762
        rin = 1.4626008660993084
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
    elseif j==2
        rcut = 5.126585004741973
        rin = 1.0407882529057981
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
    elseif j==3 
        rcut = 5.1566588199347265
        rin = 0.9
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
    elseif j==4
        rcut = 5.25
        rin = 0.9613114544981287
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)
        descriptors[2] = POD(nbody=2, pdegree=[6,8], nbasis = [11], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, pdegree=[6,8,7], nbasis = [10, 7], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
    elseif j==5        
        rcut = 5.239893145443493
        rin = 0.9815780250870151
        descriptors[1] = POD(nbody=1, pdegree=[0], nbasis = [1], rin = rin, rcut=rcut)   
        descriptors[2] = POD(nbody=2, pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
        descriptors[3] = POD(nbody=3, pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma, hybrid23=false, projectiontol=1e-10)
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

# eta = zeros(6,2)
# for j = 1:6
#     tm = readdlm("results/optpodeta" * string(j-1) *  ".txt")
#     eta[j,:] = tm[:]
# end
# show(IOContext(stdout, :compact=>false), "text/plain", eta)

