cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/InP/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["In", "P"];

# folders in the datapath 
folders = ["aa", "aIn", "aP", "Bulk",
            "EOS", "iIn", "iP", "s_aa", "s_aIn",
            "s_aP", "Shear", "s_iIn", "s_iP", "Strain", 
            "s_vIn","s_vP", "s_vv", "vP", "vv"];

n = length(folders)
weightinner = zeros(n,3)
weightinner[:,1] .= 100.0
weightinner[:,2] .= 1.0

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.8, 0.2, 0.0]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors 
transposelattice = false 

# training data 
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
                rotationmatrix, transposelattice)
end

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# bounds for cut-off radius to be optimized
rcutrange = [3.8, 5.0]

# define range for nonlinear parameters to be optimized
etarange =  reshape(rcutrange,(1,2))

# number of interpolation points for each nonlinear parameters 
N = [7]
etaspace = [etarange N]

rcut = 4.812302818   # cut-off radiuss 
wj = [1.0, 0.9590493408]
radelem = [0.5, 0.417932464]
bzeroflag = 1
chemflag = 1

# optimization parameters
eta = [rcut]; 
kappa = [0];  
optim = setoptim(lossfunc, method, eta, kappa, weightouter, etaspace)

twojmaxlist = [2, 4, 6, 8, 10]
for j = 0:4
    twojmx = twojmaxlist[j+1]

    # SNAP descriptors
    descriptors[1] = SNAPparams(species = [:In,:P], twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363, 
        elemradius = radelem, elemweight = wj, bzeroflag = bzeroflag, chemflag = chemflag)

    #j = 4 -> 4.73987269841031  0.00018358134416  0.00794100068092

    # optimize the SNAP potential 
    opteta, coeff, fmin, iter, polycoeff = optimize(traindata, traindata, descriptors, potentials, Doptions, optim)

    # compute unweighted MAE, RMSE, RSQ errors 
    energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(folders, energyerrors, "Energy Errors")
    printerrors(folders, forceerrors, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optsnapchemcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/optsnapchempolycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optsnapchemtrainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    #writedlm("results/optsnaptesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optsnapchemeta" * string(j) *  ".txt", opteta)
end

# sna2 =
#    3.989977113821321   0.019974300000000   0.115049000000000
#    4.189977113821321   0.010150000000000   0.082635000000000
#    4.389977113821321   0.008091740000000   0.069666100000000
#    4.589977113821321   0.008545250000000   0.074684700000000

# sna4 =
#    4.172630279300000   0.003534600000000   0.033495100000000
#    4.372630279300000   0.002782910000000   0.027991800000000
#    4.472630279300000   0.002743180000000   0.028104400000000
#    4.572630279300000   0.002784500000000   0.027334600000000

# sna6 =
#    4.218372560000000   0.001217980000000   0.018995000000000
#    4.418372560000000   0.001117100000000   0.017846200000000
#    4.618372560000000   0.001046870000000   0.018107100000000

# sna8 =
#    4.453297100000000   0.000437576000000   0.012689500000000
#    4.653297100000000   0.000361970000000   0.012123600000000
#    4.853297100000000   0.000366622000000   0.011799300000000

# sna10 =
#    4.853297100000000   0.000200475000000   0.008888370000000
#    4.968627000000000   0.000195148000000   0.008716930000000
#    5.068627000000000   0.000200227000000   0.008516770000000

# 4.389977113821321   0.008091740000000   0.069666100000000
# 4.572630279300000   0.002784500000000   0.027334600000000
# 4.618372560000000   0.001046870000000   0.018107100000000
# 4.853297100000000   0.000366622000000   0.011799300000000
# 5.068627000000000   0.000200227000000   0.008516770000000
