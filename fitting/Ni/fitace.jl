cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/Ni/"
dataformat = "jsonpymatgen"
fileextension = "json"
atomspecies = ["Ni"];

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

for j = 0:5    
    
    display(j)
    local acedescriptors = Array{Any}(nothing, 3)
    acedescriptors[1] = ACEpot.ACEparams(species = [:Ni], nbody=1, pdegree=1, r0=ACEpot.rnn(:Ni), rcut=5.0, rin=1.0, wL=2.0, csp=1.75)
    if j == 0
        rcut = 4.06835106522742    
        rin = 1.1667114429658787
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ni], nbody=3, pdegree=6, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        rcut = 4.387872063060331   
        rin = 1.3141741029072902
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ni], nbody=4, pdegree=8, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2
        rcut = 4.145367949010746   
        rin = 1.4205637346218714
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ni], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ni], nbody=4, pdegree=10, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        rcut = 4.080729852839264   
        rin = 1.4089442736067836
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ni], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ni], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        rcut = 3.9595942606383407  
        rin = 1.072745925528422
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ni], nbody=2, pdegree=12, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ni], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        rcut = 3.9
        rin = 1.1883286316486665
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ni], nbody=4, pdegree=13, r0=ACEpot.rnn(:Ni), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end
    
    coeff = linearfit(traindata, acedescriptors, Doptions)
    energyerrors, forceerrors = Potential.aceerroranalysis(traindata, acedescriptors, Doptions, coeff)
    energytesterrors, forcetesterrors = Potential.aceerroranalysis(testdata, acedescriptors, Doptions, coeff)

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(["train"], e1, "Energy Errors")
    printerrors(["train"], e2, "Force Errors")

    e1 = [energytesterrors[:,1] energytesterrors[:,2] 0*energytesterrors[:,1]]
    e2 = [forcetesterrors[:,1] forcetesterrors[:,2] 0*forcetesterrors[:,1]]
    printerrors(["test"], e1, "Energy Errors")
    printerrors(["test"], e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitacecoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    writedlm("results/fitacetesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
end

# eta = zeros(6,2)
# for j = 1:6
#     tm = readdlm("results/optaceeta" * string(j-1) *  ".txt")
#     eta[j,:] = tm[:]
# end
# show(IOContext(stdout, :compact=>false), "text/plain", eta)