function selectconfigurations(A1, A2, A3, Atombases, atomtype, descriptors, tol=1e-8)

E, F, G  = latticesimilaritymatrices(A1, A2, A3, Atombases, atomtype, descriptors)
#E, G  = similaritymatrices(A1, A2, A3, Atombases, atomtype, descriptors)

E = F.*G;
display(E)

if tol < 1.0
    U, W = Potential.eigenvalues(E);
    ind = findall(W .> 0.0);
    W = W[ind]
    # K = Potential.numbasis(W, tol);
    W = W/W[1];
    K = 1;
    for i = 2:length(W)
        if W[i] < tol
            K = i;
            break;
        end
    end
    display(K)
    K = min(K, 2000);
else
    K = tol
end

ind = minmaxsampling(E, K);

if size(A1,2) == 1
    lattices = [A1; A2; A3]
else
    lattices = [A1[:,ind]; A2[:,ind]; A3[:,ind]]
end

atombases =  Atombases[:,:,ind]

return lattices, atombases, ind, E

end

function lengthconfigurations(lat, ranges = nothing, dims = nothing, descriptors = nothing, tol=1e-8)

if ranges === nothing     
    amin = min(lat.a - 1.25, 1.0)
    amax = min(lat.a + 2.25, 8.0)
    bmin = min(lat.b - 1.25, 1.0)
    bmax = min(lat.b + 2.25, 8.0)
    cmin = min(lat.c - 1.25, 1.0)
    cmax = min(lat.c + 2.25, 8.0)
    ranges = [amin amax; bmin bmax; cmin cmax];        
end

if dims === nothing 
    dims = [21, 21, 21];       
end
dims = dims .- 1

if descriptors === nothing
    descriptors = PODdescriptors();
end

A1, A2, A3, Atombases = lengthgroup(lat, ranges, dims);

lattices, atombases, ind = selectconfigurations(A1, A2, A3, Atombases, lat.atomtype, descriptors, tol)

return lattices, atombases, ind, A1, A2, A3, Atombases

end

function angleconfigurations(lat, ranges = nothing, dims = nothing, descriptors = nothing, tol=1e-8)

if ranges === nothing     
    ranges = [70 110; 70 110; 70 110]
end
ranges = ranges*pi/180.0

if dims === nothing 
    dims = [21, 21, 21];       
end
dims = dims .- 1

if descriptors === nothing
    descriptors = PODdescriptors();
end

A1, A2, A3, Atombases = anglegroup(lat, ranges, dims);

lattices, atombases, ind = selectconfigurations(A1, A2, A3, Atombases, lat.atomtype, descriptors, tol)

return lattices, atombases, ind, A1, A2, A3, Atombases

end

function latticeconfigurations(lat, ranges = nothing, dims = nothing, descriptors = nothing, tol=1e-8)

if ranges === nothing     
    amin = min(lat.a - 1.0, 1.0)
    amax = min(lat.a + 1.5, 8.0)
    bmin = min(lat.b - 1.0, 1.0)
    bmax = min(lat.b + 1.5, 8.0)
    cmin = min(lat.c - 1.0, 1.0)
    cmax = min(lat.c + 1.5, 8.0)
    range1 = [amin amax; bmin bmax; cmin cmax];        
    range2 = [75 105; 75 105; 75 105];
    ranges = [range1; range2];
end
ranges[4:6,:] = ranges[4:6,:]*pi/180.0

if dims === nothing 
    dims = [6, 6, 6, 4, 4, 4];       
end
dims = dims .- 1

if descriptors === nothing
    descriptors = PODdescriptors();
end

A1, A2, A3, Atombases = latticegroup(lat, ranges, dims);

lattices, atombases, ind, E = selectconfigurations(A1, A2, A3, Atombases, lat.atomtype, descriptors, tol)

return lattices, atombases, ind, A1, A2, A3, Atombases, E

end


function atombasisconfigurations(lat, ds = nothing, dims = nothing, descriptors = nothing, tol=1e-8)

nbasis = size(lat.fraccoords,2)    
if ds === nothing         
    ds = 0.1*ones(Float64,nbasis);     
end

if dims === nothing 
    dims = 2*ones(Int64,nbasis);      
end
dims = dims .- 1

if descriptors === nothing
    descriptors = PODdescriptors();
end

A1, A2, A3, Atombases = atombasisgroup(lat, ds, dims);

lattices, atombases, ind = selectconfigurations(A1[:], A2[:], A3[:], Atombases, lat.atomtype, descriptors, tol)

return lattices, atombases, ind, A1, A2, A3, Atombases

end

function atomdisplacementconfigurations(lat, ds = nothing, N = nothing, descriptors = nothing, tol=1e-8)

nbasis = size(lat.fraccoords,2)    
if ds === nothing         
    ds = 0.1*ones(Float64,nbasis);     
end

if N === nothing
    N = 10000;        
end        

if descriptors === nothing
    descriptors = PODdescriptors();
end

A1, A2, A3, Atombases = atomdisplacementgroup(lat, ds, N);

lattices, atombases, ind = selectconfigurations(A1[:], A2[:], A3[:], Atombases, lat.atomtype, descriptors, tol)

return lattices, atombases, ind, A1, A2, A3, Atombases

end
