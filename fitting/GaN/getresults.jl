using DelimitedFiles

podeta = zeros(5,1)
podtrain = zeros(5,2)
podtest = zeros(5,2)
podnum = zeros(5)
for j = 1:5
    podcoeff = readdlm("results/fitpodcoeff" * string(j-1) *  ".txt")    
    podnum[j] = length(podcoeff)

    tm = readdlm("results/fitpodtrainerror" * string(j-1) *  ".txt")    
    podtrain[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/fitpodtesterror" * string(j-1) *  ".txt")    
    podtest[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optpodeta" * string(j-1) *  ".txt")
    podeta[j,:] = tm[:]
end
pod = [podeta podtrain podtest]
writedlm("podresults.txt", pod)

aceeta = zeros(5,1)
acetrain = zeros(5,2)
acetest = zeros(5,2)
acenum = zeros(5)
for j = 1:5
    acecoeff = readdlm("results/fitacecoeff" * string(j-1) *  ".txt")    
    acenum[j] = length(acecoeff)

    tm = readdlm("results/fitacetrainerror" * string(j-1) *  ".txt")    
    acetrain[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/fitacetesterror" * string(j-1) *  ".txt")    
    acetest[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optaceeta" * string(j-1) *  ".txt")
    aceeta[j,:] = tm[:]
end
ace = [aceeta acetrain acetest]
writedlm("aceresults.txt", ace)

pod23eta = zeros(5,1)
pod23train = zeros(5,2)
pod23test = zeros(5,2)
pod23num = zeros(5)
for j = 1:5
    pod23coeff = readdlm("results/fitpod23coeff" * string(j-1) *  ".txt")    
    pod23num[j] = length(pod23coeff)

    tm = readdlm("results/fitpod23trainerror" * string(j-1) *  ".txt")    
    pod23train[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/fitpod23testerror" * string(j-1) *  ".txt")    
    pod23test[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optpod23eta" * string(j-1) *  ".txt")
    pod23eta[j,:] = tm[:]
end
pod23 = [pod23eta pod23train pod23test]
writedlm("pod23results.txt", pod23)

snapeta = zeros(5)
snaptrain = zeros(5,2)
snaptest = zeros(5,2)
snapnum = zeros(5)
for j = 1:5
    snapcoeff = readdlm("results/fitsnapcoeff" * string(j-1) *  ".txt")    
    snapnum[j] = length(snapcoeff)

    tm = readdlm("results/fitsnaptrainerror" * string(j-1) *  ".txt")    
    snaptrain[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/fitsnaptesterror" * string(j-1) *  ".txt")    
    snaptest[j,:] = [tm[1,1], tm[1,3]]

    tm = readdlm("results/optsnapeta" * string(j-1) *  ".txt")        
    snapeta[j] = tm[1]
end
snap = [snapeta snaptrain snaptest]
writedlm("snapresults.txt", snap)

trainerror = [acetrain snaptrain podtrain pod23train]
testerror = [acetest snaptest podtest pod23test]

traine = trainerror[:,[1,3,5,7]] 
trainf = trainerror[:,[2,4,6,8]] 
teste = testerror[:,[1,3,5,7]] 
testf = testerror[:,[2,4,6,8]] 
num = [acenum snapnum podnum pod23num]

display(trainerror)
display(testerror)

writedlm("numdesc" *  ".txt", num)
writedlm("trainerror" *  ".txt", [traine trainf])
writedlm("testerror" *  ".txt", [teste testf])


