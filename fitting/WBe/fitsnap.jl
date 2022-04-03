cdir = pwd(); ii = findlast("MLP", cdir); MLPpath = cdir[1:ii[end]] * "/";    
include(MLPpath * "src/setup.jl");

using DelimitedFiles

# path to the database 
datapath = "/Users/ngoccuongnguyen/GitHub/FitSNAP/examples/WBe_PRB2019/JSON/"
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

weightinner=[41.9393253217912  68.9623173495202 0.0e-12
            41.9393253217912  68.9623173495202 0.0e-12
            41.9393253217912  68.9623173495202 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            11.6129704286138  20.1327190878775 0.0e-12
            436.46713269739    25.9225848378898 0.0e-12
            436.46713269739    25.9225848378898 0.0e-12
            178.350842653746  793.596683274656 0.0e-12
            178.350842653746 1553.7039315755 0.0e-12
            178.350842653746  793.596683274656 0.0e-12
            178.350842653746 1553.7039315755 0.0e-12
            178.350842653746  793.596683274656 0.0e-12
            178.350842653746 1553.7039315755 0.0e-12
            338.213923534691 1553.7039315755 0.0e-12
            338.213923534691 1553.7039315755 0.0e-12
            338.213923534691 1553.7039315755 0.0e-12
            11.6129704286138 793.596683274656 0.0e-12
            0.0353125583595375 25.9225848378898 0.0e-12
            44.162042355268    0.00227087929058201 0.0e-12
            44.162042355268    0.00227087929058201 0.0e-12
            0.553592561007502 0.177919040344643 0.0e-12
            211.996117831303    0.00105293249458209 0.0e-12
            132.879256296705    1.10046761515816 0.0e-12
            55.7715970380278   0.78584537072664 0.0e-12
            0.003001934       0.1396702165     0.0e-12
            0.003001934       0.1396702165     0.0e-12
            0.003001934       0.1396702165     0.0e-12
            0.003001934       0.1396702165     0.0e-12
            0.0276327087      2.1358562516     0.0e-12
            0.0001458609      0.085046817      0.0e-12
            0.0269003754      0.2379420856     0.0e-12
            0.0009078018     15.9805795873     0.0e-12
            0.0002260273      3.9286571351     0.0e-12
            0.0004723469     30.7888852447     0.0e-12
            0.0526579837      0.4206962133     0.0e-12
            0.0025503766     10.5847278636     0.0e-12
            531.0618507847    0.8123933908     0.0e-12
            3.9699030105      0.5623621081     0.0e-12]

folders = folders[25:30]            
weightinner = weightinner[25:30,:]

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

rcut = 4.812302818   # cut-off radiuss 

# single descriptors to catpure energy of isolated pure elements
#descriptors[1] = POD(nbody=1, species = [:W,:Be])

wj = [1.0, 0.9590493408]
radelem = [0.5, 0.417932464]
bzeroflag = 1
# SNAP descriptors
descriptors[1] = SNAPparams(species = [:W,:Be], twojmax = 8, rcutfac = rcut, rfac0 = 0.99363, 
    elemradius = radelem, elemweight = wj, bzeroflag = bzeroflag)

# Descriptors optional parameters
Doptions = DescriptorsOptions(pbc = [1, 1, 1], normalizeenergy=true, normalizestress=true)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforce"

# use least-square method
method = "lsq" 

# optimization parameters
optim = setoptim(lossfunc, method)

# linear fit to compute SNAP coefficients
coeff, ce, cf, cef, cefs, emae, fmae, smae,~ = linearfit(traindata, descriptors, potentials, Doptions, optim)

print("SNAP Coeffficients: "), show(stdout, "text/plain", coeff)

# # compute unweighted MAE, RMSE, RSQ erroxrs 
energyerrors, forceerrors, stresserrors = validate(traindata, descriptors, potentials, Doptions, optim, coeff)    

printerrors(folders, energyerrors, "Energy Errors")
printerrors(folders, forceerrors, "Force Errors")
printerrors(folders, stresserrors, "Stress Errors")

err = [energyerrors forceerrors stresserrors]
using DelimitedFiles
Preprocessing.mkfolder("results")
writedlm("results/fitsnapcoeff.txt", coeff)
writedlm("results/fitsnaperror.txt", err)

# ----------------------------------------------------------------------------------------
# Energy Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.07044207657701    |    0.21603314997155    |    0.99812010420952   |
# Displaced_A15 |   0.00096552143837    |    0.00111604188449    |    0.54713143273288   |
# Displaced_BCC |   0.00265790526799    |    0.00282445652492    |    0.99558694447092   |
# Displaced_FCC |   0.00062281398597    |    0.00070621825523    |    0.97482658887872   |
# Elastic_BCC   |   0.01609666719843    |    0.01611409810408    |    0.94536142997231   |
# Elastic_FCC   |   0.00236015882463    |    0.00283113393722    |    0.99731108251481   |
# GSF_110       |   0.00273706745771    |    0.00327819622044    |    0.96961211605756   |
# GSF_112       |   0.00493895723907    |    0.00543788060718    |    0.95943391379761   |
# Liquid        |   0.00202413959425    |    0.00220367987714    |    0.99954957466195   |
# Surface       |   0.00986241413964    |    0.01286642045895    |    0.99206605756607   |
# Volume_A15    |   0.15845142823669    |    0.20402247466747    |    0.99923461069253   |
# Volume_BCC    |   0.24054293169245    |    0.34857690690444    |    0.99932441035003   |
# Volume_FCC    |   0.43992415339314    |    0.65036051729788    |    0.99528215114935   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Force Errors  |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.07382175918033    |    0.13473293245264    |    0.98985401189508   |
# Displaced_A15 |   0.10121880694163    |    0.13144958315829    |    0.94699246994167   |
# Displaced_BCC |   0.11558686617518    |    0.14717456426298    |    0.99060566246513   |
# Displaced_FCC |   0.04485641678835    |    0.05799790562769    |    0.98730139228410   |
# Elastic_BCC   |   0.05124398995294    |    0.06281587640833    |    0.62885074654137   |
# Elastic_FCC   |   0.03767913962690    |    0.04833932543656    |    0.79876767281188   |
# GSF_110       |   0.02331765506779    |    0.04143778578021    |    0.99935586576560   |
# GSF_112       |   0.06227573920557    |    0.09595691088125    |    0.99759093141190   |
# Liquid        |   0.29368616066726    |    0.38087413462762    |    0.96762323524226   |
# Surface       |   0.04715202191799    |    0.10329820204475    |    0.99705018659741   |
# Volume_A15    |   3.46844637110477    |    8.17503109246182    |    0.47346230209265   |
# Volume_BCC    |   2.26883249123272    |    4.74575523467808    |    -2.5290229329517   |
# Volume_FCC    |   1.49145134984300    |    3.37372286466239    |    0.07729394577249   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Stress Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   35971.9793128237    |    195191.138554860    |    0.98421362798808   |
# Displaced_A15 |   22461.9928376593    |    31109.3032729486    |    0.99670734388147   |
# Displaced_BCC |   3156.48504229259    |    4399.91212860888    |    0.99994088314579   |
# Displaced_FCC |   11193.8066343422    |    15619.4228192125    |    0.99912377097448   |
# Elastic_BCC   |   280.397837618711    |    396.387323061027    |    0.99999951872136   |
# Elastic_FCC   |   21366.3276195083    |    29271.1138024509    |    0.99696441213386   |
# GSF_110       |   1582.59625565322    |    2337.94276204977    |    0.99993757303705   |
# GSF_112       |   2156.25006603071    |    3080.77851292823    |    0.99987925470715   |
# Liquid        |   39240.6700104592    |    54157.2641612818    |    0.98826937704021   |
# Surface       |   2273.06425492121    |    4022.91528154733    |    0.99980384280603   |
# Volume_A15    |   88790.3570675373    |    298468.446977544    |    0.96695967433505   |
# Volume_BCC    |   152032.930937683    |    452922.949353301    |    0.99110567462576   |
# Volume_FCC    |   144824.356245268    |    466410.633437398    |    0.96121748046692   |
# ----------------------------------------------------------------------------------------

