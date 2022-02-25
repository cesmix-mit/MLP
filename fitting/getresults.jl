using DelimitedFiles

# Li, Ni, Si, Cu, Mo, Ge, Ta 
folder = "Ta/"

podeta = zeros(6,2)
podtrain = zeros(6,2)
podtest = zeros(6,2)
podnum = zeros(6)
for j = 1:6
    podcoeff = readdlm(folder * "results/optpodcoeff" * string(j-1) *  ".txt")    
    podnum[j] = length(podcoeff)

    tm = readdlm(folder * "results/optpodtrainerror" * string(j-1) *  ".txt")    
    podtrain[j,:] = [tm[1,1], tm[1,4]]

    # tm = readdlm(folder * "results/optpodtesterror" * string(j-1) *  ".txt")    
    # podtest[j,:] = [tm[1,1], tm[1,4]]

    tm = readdlm(folder * "results/optpodeta" * string(j-1) *  ".txt")
    podeta[j,:] = tm[:]
end

aceeta = zeros(6,2)
acetrain = zeros(6,2)
acetest = zeros(6,2)
acenum = zeros(6)
for j = 1:6
    acecoeff = readdlm(folder * "results/optacecoeff" * string(j-1) *  ".txt")    
    acenum[j] = length(acecoeff)

    tm = readdlm(folder * "results/optacetrainerror" * string(j-1) *  ".txt")    
    acetrain[j,:] = [tm[1,1], tm[1,4]]

    # tm = readdlm(folder * "results/optacetesterror" * string(j-1) *  ".txt")    
    # acetest[j,:] = [tm[1,1], tm[1,4]]

    tm = readdlm(folder * "results/optaceeta" * string(j-1) *  ".txt")
    aceeta[j,:] = tm[:]
end

snapeta = zeros(6)
snaptrain = zeros(6,2)
snaptest = zeros(6,2)
snapnum = zeros(6)
for j = 1:6
    snapcoeff = readdlm(folder * "results/optsnapcoeff" * string(j-1) *  ".txt")    
    snapnum[j] = length(snapcoeff)

    tm = readdlm(folder * "results/optsnaptrainerror" * string(j-1) *  ".txt")    
    snaptrain[j,:] = [tm[1,1], tm[1,4]]

    # tm = readdlm(folder * "results/optsnaptesterror" * string(j-1) *  ".txt")    
    # snaptest[j,:] = [tm[1,1], tm[1,4]]

    tm = readdlm(folder * "results/optsnapeta" * string(j-1) *  ".txt")        
    snapeta[j] = tm[1]
end

trainerror = [acetrain snaptrain podtrain]
testerror = [acetest snaptest podtest]

traine = trainerror[:,[1,3,5]] 
trainf = trainerror[:,[2,4,6]] 
teste = testerror[:,[1,3,5]] 
testf = testerror[:,[2,4,6]] 
num = [acenum snapnum podnum]

display(trainerror)
display(testerror)

writedlm(folder * "numdesc" *  ".txt", num)
writedlm(folder * "trainerror" *  ".txt", [traine trainf])
writedlm(folder * "testerror" *  ".txt", [teste testf])


