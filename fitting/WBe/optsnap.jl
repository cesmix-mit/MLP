cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "../../data/WBe/JSON/"
dataformat = "json"
fileextension = "json"
atomspecies = ["W", "Be"];

# folders in the datapath 
folders = ["001FreeSurf", "010FreeSurf", "100FreeSurf", "Defect_BCrowd",
            "Defect_BOct", "Defect_BSplit", "Defect_BTetra", "Defect_Crowd", "Defect_Oct",
            "Defect_Tet", "Defect_Vacancy", "DFTMD_1000K", "DFTMD_300K", "Elast_BCC_Shear", 
            "Elast_BCC_Vol","Elast_FCC_Shear", "Elast_FCC_Vol", "Elast_HCP_Shear", "Elast_HCP_Vol",
            "EOS_BCC","EOS_FCC", "EOS_HCP", "StackFaults","Liquids","DFT_MD_1000K", "DFT_MD_300K", 
            "ElasticDeform_Shear","ElasticDeform_Vol", "EOS", "WSurface_BeAdhesion", 
            "BCC_ForceLib_W110", "BCC_ForceLib_W111", "BCC_ForceLib_WOct", "BCC_ForceLib_WTet",
            "dislocation_quadrupole", "Disordered_Struc", "Divacancy", "EOS_Data", "gamma_surface",
            "gamma_surface_vacancy","md_bulk","slice_sample","surface","vacancy"];

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
chemflag = 0

# optimization parameters
eta = [rcut]; 
kappa = [0];  
optim = setoptim(lossfunc, method, eta, kappa, weightouter, etaspace)

twojmaxlist = [2, 4, 6, 8, 10]
for j = 1:1
    twojmx = twojmaxlist[j+1]

    # SNAP descriptors
    descriptors[1] = SNAPparams(species = [:W,:Be], twojmax = twojmx, rcutfac = rcut, rfac0 = 0.99363, 
        elemradius = radelem, elemweight = wj, bzeroflag = bzeroflag, chemflag=chemflag)

    # optimize the SNAP potential 
    opteta, coeff, fmin, iter, polycoeff = optimize(traindata, traindata, descriptors, potentials, Doptions, optim)

    # compute unweighted MAE, RMSE, RSQ errors 
    energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(folders, energyerrors, "Energy Errors")
    printerrors(folders, forceerrors, "Force Errors")

    Preprocessing.mkfolder("results")
    writedlm("results/optsnap2coeff" * string(j) *  ".txt", coeff)
    writedlm("results/optsnap2polycoeff" * string(j) *  ".txt", polycoeff)
    writedlm("results/optsnap2trainerror" * string(j) *  ".txt", [energyerrors forceerrors])    
    #writedlm("results/optsnaptesterror" * string(j) *  ".txt", [energytesterrors forcetesterrors])    
    writedlm("results/optsnap2eta" * string(j) *  ".txt", opteta)
end

# snap2 =
#     3.8000    0.1246    0.5280
#     4.1266    0.1613    0.6590

# snap4 =
#     3.8000    0.0901    0.2359
#     3.9609    0.0822    0.2110
#     4.2139    0.0780    0.1999
#     4.4028    0.0801    0.2249

# snap6 =
#     4.4028    0.0702    0.1961
#     4.6193    0.0706    0.2432
#     4.2193    0.0713    0.1779

# snap8 =
#     4.1938    0.0696    0.1839
#     4.3938    0.0687    0.1992
#     4.2938    0.0692    0.1855


