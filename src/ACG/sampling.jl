function minmaxsampling(S, k)
# S: similarity matrix 
    
n = size(S,1);

amin = zeros(n,1);
jmin = zeros(n,1);
b = zeros(k,1);
ind = zeros(Int32,k,1);

for i = 1:n
    amin[i], jmin[i] = findmin(S[:,i]);    
end
tm, imin = findmin(amin);

ind[1] = imin[1];
ind[2] = jmin[imin];

for p = 3:k
    indp = ind[1:(p-1)];
    m = p - 1;
    for i = 1:n
        for j = 1:m
            b[j] = S[i,indp[j]];            
        end
        amin[i] = maximum(b[1:m]);        
    end
    tm, imin = findmin(amin);
    ind[p] = imin[1];
end

return ind

end

function maxminsampling(D, k)
# D: dissimilarity matrix 

n = size(D,1);

amax = zeros(n,1);
jmax = zeros(n,1);
b = zeros(k,1);
ind = zeros(Int32,k,1);

for i = 1:n
    amax[i], jmax[i] = findmax(D[:,i]);    
end
tm, imax = findmax(amax);

ind[1] = imax;
ind[2] = jmax[imax];

for p = 3:k
    indp = ind[1:(p-1)];
    m = p - 1;
    for i = 1:n
        for j = 1:m
            b[j] = D[i,indp[j]];            
        end
        amax[i] = minimum(b[1:m]);        
    end
    tm, imax = findmax(amax);
    ind[p] = imax;
end

return ind

end

