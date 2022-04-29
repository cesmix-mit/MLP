cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "installation/setup.jl");

# path to the database 
datapath = "../../data/Ta/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["Ta"];

# folders in the datapath 
folders = ["Displaced_A15", "Displaced_BCC", "Displaced_FCC", "Elastic_BCC",
            "Elastic_FCC", "GSF_110", "GSF_112", "Liquid", "Surface",
            "Volume_A15", "Volume_BCC", "Volume_FCC"];

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

weightinner = nothing

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = true 

# training data 
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner, translationvector, 
                rotationmatrix, transposelattice)
end

rcut = 4.714999753078612   # cut-off radiuss 

# single descriptors to catpure energy of isolated pure elements
descriptors[1] = PotentialStruct("single", "nonbonded", [1], "single", [0.0], rcut);
include(descriptors[1].potentialfunction * ".jl");  # include the descriptors file

# SNAP descriptors
descriptors[2] = initsnap(SNAPparams(species = [:Ta], twojmax = 8, rcutfac = rcut, rfac0 = 0.99363))

# ZBL reference potential
potentials[1] = addpotential("pair", "nonbonded", [1, 1], "zbl", [4.0, 4.8, 73], 4.8);
include(potentials[1].potentialfunction * ".jl");  # include the potential file

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=false, normalizestress=false)

# activation function 
# f(x) = (1.0)./(1.0 .+ exp(-x)) - 0.5

# # Derivative of activation function 
# df(x) = exp(-x)./((exp(-x) + 1).^2) # f(x).*(1-f(x)) 

# activation function 
f(x) = (2.0)./(1.0 .+ exp(-0.25*x)) - 1.0

# Derivative of activation function 
df(x) = 0.5*exp(-0.25*x)./((exp(-0.25*x) + 1).^2) 

# f = (1.0)./(1.0 + exp(-a*x))
# df = (a*exp(-a*x))/(exp(-a*x) + 1)^2

# Set model parameters 
η = [1e-3, 1e-3, 1e-3]
modelOpts = DML.ModelOptions(modeltype = "dnn", lossfunc = DML.lossmse, f = f, df = df, η = η, 
                 epochs = 5000, num_layers = 3, batchsize = 16, β = [0.2, 0.8])

#  train DNN model 
model, data = DML.TrainModel(traindata, descriptors, Doptions, modelOpts, potentials)                  

# calculate errors 
maeenergy, maeforce = DML.MAE(model, data, modelOpts)

# # save the model 
using JLD2
Preprocessing.mkfolder("results")
jldsave("results/dnnsnapmodel.jld2"; model, modelOpts, data)

# tm = load("dnnsnapmodel.jld2")
# model = tm["model"]

err = [maeenergy maeforce]
using DelimitedFiles
writedlm("results/dnnsnaperror.txt", err)

# Energy MAE Error: 0.09385775263014601 
# Force MAE Error: 0.06883991876823732 
# "Loss value: 5.2287321373433695"

# Energy MAE Error: 0.013561782233723584 
# Force MAE Error: 0.05541929867207133 
# "Loss value: 0.6491575931224266"

# Energy MAE Error: 0.011437660253318304 
# Force MAE Error: 0.05337999428665492 
# "Loss value: 0.5937190723052249"

