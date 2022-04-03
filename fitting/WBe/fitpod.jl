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

# folders = folders[29:30]
# weightinner = weightinner[29:30,:]

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

# optimization parameters
optim = setoptim(lossfunc, method)

# cut-off radius
rcut = 4.691630445488325

# inner radius 
rin = 0.9

# Bessel scaling parameters
gamma = [0.0, 2, 4]
for j = 0:0

    display(j)

    # POD Descriptors
    if j == 0
        descriptors[1] = POD(nbody=2, species = [:W,:Be], pdegree=[2,4], nbasis = [2], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:W,:Be], pdegree=[2,2,2], nbasis = [1, 2], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j == 1
        descriptors[1] = POD(nbody=2, species = [:W,:Be], pdegree=[2,6], nbasis = [3], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:W,:Be], pdegree=[2,4,3], nbasis = [3, 3], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==2
        descriptors[1] = POD(nbody=2, species = [:W,:Be], pdegree=[3,6], nbasis = [6], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:W,:Be], pdegree=[3,6,4], nbasis = [5, 4], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==3 
        descriptors[1] = POD(nbody=2, species = [:W,:Be], pdegree=[4,6], nbasis = [8], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:W,:Be], pdegree=[4,6,5], nbasis = [8, 5], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==4
        descriptors[1] = POD(nbody=2, species = [:W,:Be], pdegree=[6,8], nbasis = [14], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:W,:Be], pdegree=[6,8,7], nbasis = [11, 6], rin = rin, rcut=rcut, gamma0 = gamma)
    elseif j==5           
        descriptors[1] = POD(nbody=2, species = [:W,:Be], pdegree=[6,12], nbasis = [10], rin = rin, rcut=rcut, gamma0 = gamma)
        descriptors[2] = POD(nbody=3, species = [:W,:Be], pdegree=[6,12,10], nbasis = [12,10], rin = rin, rcut=rcut, gamma0 = gamma)
    end

    # linear fit to compute SNAP coefficients
    coeff, ce, cf, cef, cefs, emae, fmae, smae,~ = linearfit(traindata, descriptors, potentials, Doptions, optim)

    print("POD Coeffficients: "), show(stdout, "text/plain", coeff)

    # # compute unweighted MAE, RMSE, RSQ erroxrs 
    energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

    printerrors(folders, energyerrors, "Energy Errors")
    printerrors(folders, forceerrors, "Force Errors")
    printerrors(folders, stresserrors, "Stress Errors")

    err = [energyerrors forceerrors stresserrors]
    using DelimitedFiles
    Preprocessing.mkfolder("results")
    writedlm("results/fitpodcoeff" * string(j) *  ".txt", coeff)
    writedlm("results/fitpoderror" * string(j) *  ".txt", err)
end

