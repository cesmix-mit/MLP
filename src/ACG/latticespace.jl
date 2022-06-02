using Optimization 

function lengthspace(ranges=nothing, dims=nothing)

    if ranges === nothing
        ranges = [1.5 5; 1.5 5; 1.5 5];        
    end    
    if dims === nothing
        dims = [25, 25, 25].-1;        
    end    
    
    return Optimization.tensornodes(ranges[:,1], ranges[:,2], dims)
    
end

function anglespace(ranges=nothing, dims=nothing)
    
    if ranges === nothing
        ranges = [45 135; 45 135; 45 135]*pi/180.0;        
    end    
    if dims === nothing
        dims = [25, 25, 25].-1;        
    end    
    
    return Optimization.tensornodes(ranges[:,1], ranges[:,2], dims)
    
end

function latticespace(ranges=nothing, dims=nothing)
    
    if ranges === nothing
        ranges = [1.7 4.9; 1.7 4.9; 1.7 4.9; 60 120; 60 120; 60 120];        
    end    
    if dims === nothing
        dims = [20, 20, 20, 4, 4, 4].-1;        
    end    
    
    return Optimization.tensornodes(ranges[:,1], ranges[:,2], dims)    
end

function atombasisspace(fraccoords, ds=nothing, dims=nothing)
    
    nbasis = size(fraccoords,2)
    if ds === nothing
        ds = 0.1*ones(Float64,nbasis);     
    end    
    if dims === nothing
        dims = 4*ones(Int64,nbasis);        
    end    
    
    atomspace = Array{Any}(nothing, nbasis)
    for n = 1:nbasis  
        if dims[n] > 1
            s = fraccoords[:,n]
            smin = s .- ds[n]
            smax = s .+ ds[n]
            atomspace[n] = Optimization.tensornodes(smin, smax, [dims[n], dims[n], dims[n]])
            M = length(atomspace[n][:,1])
            #display(size(atomspace[n]))
            for d = 1:3
                for j = 1:M
                    if atomspace[n][j,d] < 0.0
                        atomspace[n][j,d] = atomspace[n][j,d] + 1.0
                    elseif atomspace[n][j,d] > 1.0
                        atomspace[n][j,d] = atomspace[n][j,d] - 1.0   
                    end     
                end
            end
        else        
            atomspace[n] = zeros(1,3)
            atomspace[n][1,1] = fraccoords[1,n]
            atomspace[n][1,2] = fraccoords[2,n]
            atomspace[n][1,3] = fraccoords[3,n]
        end    
    end
    
    return atomspace
    
end
    

function atombasis(fraccoords, ds=nothing, N=nothing)

    
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
        AtomBasis[:,:,n] = w
    end

    return AtomBasis
end

