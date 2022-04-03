using DelimitedFiles

podeta = zeros(6,2)
podtrain = zeros(6,2)
podtest = zeros(6,2)
podnum = zeros(6)
podtrainerror = Array{Any}(nothing, 6)
for j = 1:6
    podcoeff = readdlm("results/fitpodcoeff" * string(j-1) *  ".txt")    
    podnum[j] = length(podcoeff)

    tm = readdlm("results/fitpodtrainerror" * string(j-1) *  ".txt")    
    podtrain[j,:] = [tm[1,1], tm[1,3]]
    podtrainerror[j] = tm[1:end,[1, 3]]

    # tm = readdlm("results/fitpodtesterror" * string(j-1) *  ".txt")    
    # podtest[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optpodeta" * string(j-1) *  ".txt")
    podeta[j,:] = tm[:]
end

aceeta = zeros(6,2)
acetrain = zeros(6,2)
acetest = zeros(6,2)
acenum = zeros(6)
acetrainerror = Array{Any}(nothing, 6)
for j = 1:6
    acecoeff = readdlm("results/fitacecoeff" * string(j-1) *  ".txt")    
    acenum[j] = length(acecoeff)

    tm = readdlm("results/fitacetrainerror" * string(j-1) *  ".txt")    
    acetrain[j,:] = [tm[1,1], tm[1,3]]
    acetrainerror[j] = tm[1:end,[1, 3]]

    # tm = readdlm("results/fitacetesterror" * string(j-1) *  ".txt")    
    # acetest[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optaceeta" * string(j-1) *  ".txt")
    aceeta[j,:] = tm[:]
end

snapeta = zeros(6)
snaptrain = zeros(6,2)
snaptest = zeros(6,2)
snapnum = zeros(6)
snaptrainerror = Array{Any}(nothing, 6)
for j = 1:6
    snapcoeff = readdlm("results/fitsnapcoeff" * string(j-1) *  ".txt")    
    snapnum[j] = length(snapcoeff)

    tm = readdlm("results/fitsnaptrainerror" * string(j-1) *  ".txt")    
    snaptrain[j,:] = [tm[1,1], tm[1,3]]
    snaptrainerror[j] = tm[1:end,[1, 3]]

    # tm = readdlm("results/fitsnaptesterror" * string(j-1) *  ".txt")    
    # snaptest[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optsnapeta" * string(j-1) *  ".txt")        
    snapeta[j] = tm[1]
end

trainerror = [acetrain snaptrain podtrain]
testerror = [acetest snaptest podtest]

traine = trainerror[:,[1,3,5]] 
trainf = trainerror[:,[2,4,6]] 
teste = testerror[:,[1,3,5]] 
testf = testerror[:,[2,4,6]] 
num = [acenum snapnum podnum]

e2 = 1e3*[acetrainerror[2] snaptrainerror[2] podtrainerror[2]] 
e4 = 1e3*[acetrainerror[4] snaptrainerror[4] podtrainerror[4]] 
e6 = 1e3*[acetrainerror[6] snaptrainerror[6] podtrainerror[6]] 

display(trainerror)
display(testerror)

writedlm("numdesc" *  ".txt", num)
writedlm("trainerror" *  ".txt", [traine trainf])
writedlm("testerror" *  ".txt", [teste testf])


