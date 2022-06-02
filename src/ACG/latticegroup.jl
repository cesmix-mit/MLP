function lengthgroup(a, b, c, alpha, beta, gamma, fraccoords, ranges=nothing, dims=nothing)
    
    cosalpha = cos(alpha);
    cosbeta = cos(beta);
    cosgamma = cos(gamma);    
    space = lengthspace(ranges, dims)
    lx, ly, lz, xy, xz, yz = triclinicbox(space[:,1], space[:,2], space[:,3], cosalpha, cosbeta, cosgamma)
        
    return tricliniclattice(lx, ly, lz, xy, xz, yz, fraccoords)
end

function lengthgroup(cosalpha, cosbeta, cosgamma, fraccoords, ranges=nothing, dims=nothing)
    
    space = lengthspace(ranges, dims)
    lx, ly, lz, xy, xz, yz = triclinicbox(space[:,1], space[:,2], space[:,3], cosalpha, cosbeta, cosgamma)
        
    return tricliniclattice(lx, ly, lz, xy, xz, yz, fraccoords)
end

function lengthgroup(lat, ranges=nothing, dims=nothing)
    
    space = lengthspace(ranges, dims)
    lx, ly, lz, xy, xz, yz = triclinicbox(space[:,1], space[:,2], space[:,3], lat.cosalpha, lat.cosbeta, lat.cosgamma)
        
    return tricliniclattice(lx, ly, lz, xy, xz, yz, lat.fraccoords)
end

function anglegroup(a, b, c, alpha, beta, gamma, fraccoords, ranges=nothing, dims=nothing)
        
    space = anglespace(ranges, dims)
    cosalpha = cos.(space[:,1]);
    cosbeta = cos.(space[:,2]);
    cosgamma = cos.(space[:,3]);
    lx, ly, lz, xy, xz, yz = triclinicbox(a, b, c, cosalpha, cosbeta, cosgamma)
    return tricliniclattice(lx, ly, lz, xy, xz, yz, fraccoords)
end

function anglegroup(a, b, c, fraccoords, ranges=nothing, dims=nothing)
        
    space = anglespace(ranges, dims)
    cosalpha = cos.(space[:,1]);
    cosbeta = cos.(space[:,2]);
    cosgamma = cos.(space[:,3]);
    lx, ly, lz, xy, xz, yz = triclinicbox(a, b, c, cosalpha, cosbeta, cosgamma)
    return tricliniclattice(lx, ly, lz, xy, xz, yz, fraccoords)
end

function anglegroup(lat, ranges=nothing, dims=nothing)
        
    space = anglespace(ranges, dims)
    cosalpha = cos.(space[:,1]);
    cosbeta = cos.(space[:,2]);
    cosgamma = cos.(space[:,3]);        
    lx, ly, lz, xy, xz, yz = triclinicbox(lat.a, lat.b, lat.c, cosalpha, cosbeta, cosgamma)
    return tricliniclattice(lx, ly, lz, xy, xz, yz, lat.fraccoords)
end

function latticegroup(a, b, c, alpha, beta, gamma, fraccoords, ranges=nothing, dims=nothing)
        
    space = latticespace(ranges, dims)    
    cosalpha = cos.(space[:,4]);
    cosbeta = cos.(space[:,5]);
    cosgamma = cos.(space[:,6]);
    lx, ly, lz, xy, xz, yz = triclinicbox(space[:,1], space[:,2], space[:,3], cosalpha, cosbeta, cosgamma)
    return tricliniclattice(lx, ly, lz, xy, xz, yz, fraccoords)
end

function latticegroup(lat, ranges=nothing, dims=nothing)
        
    space = latticespace(ranges, dims)
    cosalpha = cos.(space[:,4]);
    cosbeta = cos.(space[:,5]);
    cosgamma = cos.(space[:,6]);
    lx, ly, lz, xy, xz, yz = triclinicbox(space[:,1], space[:,2], space[:,3], cosalpha, cosbeta, cosgamma)
    return tricliniclattice(lx, ly, lz, xy, xz, yz, lat.fraccoords)
end

function atombasisgroup(a, b, c, alpha, beta, gamma, fraccoords, ds=nothing, dims=nothing)
    
    cosalpha = cos(alpha);
    cosbeta = cos(beta);
    cosgamma = cos(gamma);    
    a1, a2, a3, lx, ly, lz, xy, xz, yz   = tricliniclatticevectors(a, b, c, cosalpha, cosbeta, cosgamma)
    
    space = atombasisspace(fraccoords, ds, dims)    
    nbasis = size(fraccoords,2)  
    Ns = ones(nbasis)
    for n = 1:nbasis
        Ns[n] = length(space[n][1])
    end

    Na = ones(nbasis,1);
    for m = 2:nbasis
        Na[m] = prod(Ns[1:(m-1)]);
    end

    N = Int64(prod(Ns))    
    AtomBasis = zeros(3,nbasis,N)
    for n = 1:N
        for m = 1:nbasis
            i = Int64(rem(ceil(n/Na[m])-1,Ns[m])+1);
            z1 = space[m][i,1]
            z2 = space[m][i,2]
            z3 = space[m][i,3]
            AtomBasis[1,m,n] = lx*z1 + xy*z2 + xz*z3 
            AtomBasis[2,m,n] = ly*z2 + yz*z3      
            AtomBasis[3,m,n] = lz*z3             
        end    
    end
    
    return a1, a2, a3, AtomBasis  
end

function atombasisgroup(lat, ds=nothing, dims=nothing)

    a1, a2, a3, lx, ly, lz, xy, xz, yz = tricliniclatticevectors(lat.a, lat.b, lat.c, lat.cosalpha, lat.cosbeta, lat.cosgamma)
    
    space = atombasisspace(lat.fraccoords, ds, dims)    
    nbasis = size(lat.fraccoords,2)  
    Ns = ones(nbasis)
    for n = 1:nbasis
        Ns[n] = size(space[n],1)
    end

    Na = ones(nbasis,1);
    for m = 2:nbasis
        Na[m] = prod(Ns[1:(m-1)]);
    end

    N = Int64(prod(Ns))    
    AtomBasis = zeros(3,nbasis,N)
    for n = 1:N
        for m = 1:nbasis
            i = Int64(rem(ceil(n/Na[m])-1,Ns[m])+1);
            z1 = space[m][i,1]
            z2 = space[m][i,2]
            z3 = space[m][i,3]
            AtomBasis[1,m,n] = lx*z1 + xy*z2 + xz*z3 
            AtomBasis[2,m,n] = ly*z2 + yz*z3      
            AtomBasis[3,m,n] = lz*z3             
        end    
    end
    
    return a1, a2, a3, AtomBasis 

end


function atomdisplacementgroup(lat, ds=nothing, N=nothing)

    a1, a2, a3, lx, ly, lz, xy, xz, yz = tricliniclatticevectors(lat.a, lat.b, lat.c, lat.cosalpha, lat.cosbeta, lat.cosgamma)
    
    fraccoords = lat.fraccoords;
    nbasis = size(fraccoords,2)
    if ds === nothing
        ds = 0.1*ones(Float64,nbasis);     
    end    
    ds = [ds ds ds]';

    if N === nothing
        N = 10000;        
    end        

    AtomBasis = zeros(3,nbasis,N)
    for n = 1:N
        v = 2*rand(3, nbasis) .- 1.0
        w = fraccoords - v .* ds
        for d = 1:3
            for j = 1:nbasis
                if w[d,j] < 0.0
                    w[d,j] = w[d,j] + 1.0
                elseif w[d,j] > 1.0
                    w[d,j] = w[d,j] - 1.0
                end     
            end
        end
        for j = 1:nbasis
            z1 = w[1,j]
            z2 = w[2,j]
            z3 = w[3,j]
            AtomBasis[1,j,n] = lx*z1 + xy*z2 + xz*z3 
            AtomBasis[2,j,n] = ly*z2 + yz*z3      
            AtomBasis[3,j,n] = lz*z3             
        end
    end
    
    return a1, a2, a3, AtomBasis 
end
