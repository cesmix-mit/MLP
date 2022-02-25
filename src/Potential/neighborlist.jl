
function neighborlist(r, rcut)

dim, N = size(r);

r1 = reshape(r,(dim, N, 1));
r2 = reshape(r,(dim, 1, N));
r12 = r1 .- r2;
d = reshape(sqrt.(sum(r12.*r12,dims=1)), (N,N)); # distances between atoms

A = (d .<= rcut);  # adjacency marix
index = 1:(N+1):N*N; # indices of diagonal entries
A[index] .= 0.0;  # set diagonal of A to zero

# pairs of atoms (i,j) within the cut-off radius rcut
ind = findall(A.==1);
ai = getindex.(ind, 2);
aj = getindex.(ind, 1);

# number of neighbors for each atom i
numneigh = accumarray(ai, ones(length(ai),1)); 

ai = Int64.(ai);
aj = Int64.(aj);
numneigh = Int64.(numneigh);

return ai, aj, numneigh

end

function neighborlist(r, rcut, nx)

    rcutsq = rcut*rcut
    dim, N = size(r);
    numneigh = zeros(Int64,nx)
    ai = zeros(Int64,nx*N)
    aj = zeros(Int64,nx*N)
    k = 0
    for i = 1:nx
        ri = r[:,i]
        inc = 0
        for j = 1:N
            rj = r[:,j]
            rij = ri - rj
            rijsq = rij[1]*rij[1] + rij[2]*rij[2] + rij[3]*rij[3]
            if  (rijsq > 1e-12) & (rijsq <= rcutsq)    
                inc += 1           
                k += 1           
                ai[k] = i
                aj[k] = j                                                      
            end
        end
        numneigh[i] = inc 
    end

    ai = ai[1:k]
    aj = aj[1:k]

    return ai, aj, numneigh

end



