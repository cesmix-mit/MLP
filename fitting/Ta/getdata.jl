cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

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
config = Array{Any}(nothing, length(folders))
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
                rotationmatrix, transposelattice)

    config[i], indices = Preprocessing.readconfigdata(traindata[i])                                  
end


e = config[11].e[:]./config[11].natom[:]
v = (config[11].lattice[1,:].*config[11].lattice[5,:].*config[11].lattice[9,:])./config[11].natom[:]
bcc = [e v]
show(stdout, "text/plain", bcc)

e = config[12].e[:]./config[12].natom[:]
v = (config[12].lattice[1,:].*config[12].lattice[5,:].*config[12].lattice[9,:])./config[12].natom[:]
fcc = [e v]
show(stdout, "text/plain", fcc)

e = config[10].e[:]./config[10].natom[:]
v = (config[10].lattice[1,:].*config[10].lattice[5,:].*config[10].lattice[9,:])./config[10].natom[:]
a15 = [e v]
show(stdout, "text/plain", a15)
