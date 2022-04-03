cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/Ta/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["Ta"];

# folders in the datapath 
folders = ["Displaced_A15", "Displaced_BCC", "Displaced_FCC", "Elastic_BCC",
            "Elastic_FCC", "GSF_110", "GSF_112", "Liquid", "Surface",
            "Volume_A15", "Volume_BCC", "Volume_FCC"];

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

eta = zeros(6,2)
for j = 1:6
    tm = readdlm("results/optaceeta" * string(j-1) *  ".txt")
    eta[j,:] = tm[:]
end
# show(IOContext(stdout, :compact=>false), "text/plain", eta)

for j = 0:5    
    
    rcut = eta[j+1,1]
    rin = eta[j+1,2]
    display([j rcut rin])

    local acedescriptors = Array{Any}(nothing, 3)
    acedescriptors[1] = ACEpot.ACEparams(species = [:Ta], nbody=1, pdegree=1, r0=ACEpot.rnn(:Ta), rcut=5.0, rin=1.0, wL=2.0, csp=1.75)
    if j == 0
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=3, pdegree=6, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=2.0, csp=1.75)
    elseif j == 1
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=8, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.75, csp=1.5)
    elseif j==2
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=10, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.25, csp=1.5)
    elseif j==3 
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=2, pdegree=3, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
    elseif j==4
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=2, pdegree=12, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.5, csp=1.5)
        acedescriptors[3] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=12, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.35, csp=1.25)
    elseif j==5           
        acedescriptors[2] = ACEpot.ACEparams(species = [:Ta], nbody=4, pdegree=13, r0=ACEpot.rnn(:Ta), rcut=rcut, rin=rin, wL=1.15, csp=1.25)
    end
    
    coeff = linearfit(traindata, acedescriptors, Doptions)
    energyerrors, forceerrors, stresserrors, eerr, ferr, ferm, energies = Potential.aceerroranalysis(traindata, acedescriptors, Doptions, coeff)
    
    show(stdout, "text/plain", energies[11])

    e1 = [energyerrors[:,1] energyerrors[:,2] 0*energyerrors[:,1]]
    e2 = [forceerrors[:,1] forceerrors[:,2] 0*forceerrors[:,1]]
    printerrors(folders, e1, "Energy Errors")
    printerrors(folders, e2, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/fitacecoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitacetrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
end

