
ndgrid(v::AbstractVector) = (copy(v),)

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j - 1, snext), s) + 1]
    end
end

function ndgrid(vs::AbstractVector{T}...) where {T}
    n = length(vs)
    sz = map(length, vs)    
    out = ntuple(i -> Array{T}(undef, sz), n)
    s = 1
    for i = 1:n
        a = out[i]::Array
        v = vs[i]
        snext = s * size(a, i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

function xchenodes(N)
    # Extended Chebyshev nodes on the interval [-1, 1]

    x = 0.0;
    if N>0    
        n = N+1; 
        k = Array(1:n);       
        x = -cos.((2*k .- 1)*pi/(2*n))/cos(pi/(2*n));    
        #x = 0.5 .+ 0.5*x;
    end
    return x
end    

function tensornodes(N::Vector{Int64})
    
    dim = length(N)
    if dim==1
        x = xchenodes(N[1])
        return x
    end

    n = prod(N .+ 1)
    x1 = xchenodes(N[1])
    x2 = xchenodes(N[2])
        
    if dim==2
        xt = ndgrid(x1, x2)        
    elseif dim==3
        x3 = xchenodes(N[3])
        xt = ndgrid(x1, x2, x3)        
    elseif dim==4
        x3 = xchenodes(N[3])
        x4 = xchenodes(N[4])
        xt = ndgrid(x1, x2, x3, x4)        
    elseif dim==5
        x3 = xchenodes(N[3])
        x4 = xchenodes(N[4])
        x5 = xchenodes(N[5])
        xt = ndgrid(x1, x2, x3, x4, x5)        
    elseif dim==6
        x3 = xchenodes(N[3])
        x4 = xchenodes(N[4])
        x5 = xchenodes(N[5])
        x6 = xchenodes(N[6])
        xt = ndgrid(x1, x2, x3, x4, x5, x6)        
    elseif dim==7
        xi = Array{Any}(undef, dim)    
        for i = 1:dim
            xi[i] = xchenodes(N[i])
        end
        xt = ndgrid(xi[1], xi[2], xi[3], xi[4], xi[5], xi[6], xi[7])      
    elseif dim==8
        xi = Array{Any}(undef, dim)    
        for i = 1:dim
            xi[i] = xchenodes(N[i])
        end
        xt = ndgrid(xi[1], xi[2], xi[3], xi[4], xi[5], xi[6], xi[7], xi[8])      
    else
        error("Dimension must not exceed 8")  
    end

    x = zeros(n, dim)        
    for i = 1:dim
        x[:,i] = xt[i][:]
    end

    return x
end


function legendrepoly(x::Vector{Float64}, N::Int64)
    
    m = length(x)
    k = N + 1
    p = ones(m,k)

    if N > 0
        p[:,2] = x        
    end

    if N > 1
        for i = 3:k
            p[:,i] = ((2*(i-2)+1)/(i-1))*(x.*p[:,i-1]) - ((i-2)/(i-1))*p[:,i-2]
        end
    end

    return p
end

function legendrepolyder(x::Vector{Float64}, N::Int64)
    
    m = length(x)
    k = N + 1
    p = ones(m,k)
    dp = zeros(m, k)

    if N > 0
        p[:,2] = x      
        dp[:,2] .= 1.0  
    end

    if N > 1
        for i = 3:k
            p[:,i] = ((2*(i-2)+1)/(i-1))*(x.*p[:,i-1]) - ((i-2)/(i-1))*p[:,i-2]
            dp[:,i] = ((2*(i-2)+1)/(i-1))*(p[:,i-1] + x.*dp[:,i-1]) - ((i-2)/(i-1))*dp[:,i-2]
        end
    end

    return p, dp
end

function tensorpoly(x::Matrix{Float64}, N::Vector{Int64})
    
    n,dim = size(x)
    if dim==1
        p = legendrepoly(x[:,1], N[1])
        return p 
     end
 
    p = zeros(n, prod(N .+ 1))    
    p1 = legendrepoly(x[:,1], N[1])
    p2 = legendrepoly(x[:,2], N[2])    
    if dim==2
        for i = 1:n            
            p[i,:] = kron(p1[i,:], p2[i,:])
        end
    elseif dim==3
        p3 = legendrepoly(x[:,3], N[3])    
        for i = 1:n
            p[i,:] = kron(kron(p1[i,:], p2[i,:]), p3[i,:])
        end
    elseif dim==4
        p3 = legendrepoly(x[:,3], N[3])    
        p4 = legendrepoly(x[:,4], N[4])    
        for i = 1:n
            p[i,:] = kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:])
        end    
    elseif dim==5
        p3 = legendrepoly(x[:,3], N[3])    
        p4 = legendrepoly(x[:,4], N[4])   
        p5 = legendrepoly(x[:,5], N[5])    
        for i = 1:n
            p[i,:] = kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:])
        end          
    elseif dim==6
        p3 = legendrepoly(x[:,3], N[3])    
        p4 = legendrepoly(x[:,4], N[4])   
        p5 = legendrepoly(x[:,5], N[5])    
        p6 = legendrepoly(x[:,6], N[6])    
        for i = 1:n
            p[i,:] = kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[:,i])
        end                  
    elseif dim==7
        p3 = legendrepoly(x[:,3], N[3])    
        p4 = legendrepoly(x[:,4], N[4])   
        p5 = legendrepoly(x[:,5], N[5])    
        p6 = legendrepoly(x[:,6], N[6])    
        p7 = legendrepoly(x[:,7], N[7])    
        for i = 1:n
            p[i,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[:,i]), p7[:,i])
        end                          
    else
        error("Dimension must not exceed 7")  
    end

    return p
end

function tensorpolyder(x::Matrix{Float64}, N::Vector{Int64})
    
    n,dim = size(x)
    if dim==1
        p, dp = legendrepolyder(x[:,1], N[1])
        return p, dp 
     end

    p = zeros(n, prod(N .+ 1))
    dp = zeros(n, dim, prod(N .+ 1))
    
    p1, dpx1 = legendrepolyder(x[:,1], N[1])
    p2, dpx2 = legendrepolyder(x[:,2], N[2])    
    if dim==2
        for i = 1:n            
            p[i,:] = kron(p1[i,:], p2[i,:])
            dp[i,1,:] = kron(dpx1[i,:], p2[i,:])
            dp[i,2,:] = kron(p1[i,:], dpx2[i,:])
        end
    elseif dim==3
        p3, dpx3 = legendrepolyder(x[:,3], N[3])    
        for i = 1:n
            p[i,:] = kron(kron(p1[i,:], p2[i,:]), p3[i,:])
            dp[i,1,:] = kron(kron(dpx1[i,:], p2[i,:]), p3[i,:])
            dp[i,2,:] = kron(kron(p1[i,:], dpx2[i,:]), p3[i,:])
            dp[i,3,:] = kron(kron(p1[i,:], p2[i,:]), dpx3[i,:])
        end
    elseif dim==4
        p3, dpx3 = legendrepolyder(x[:,3], N[3])    
        p4, dpx4 = legendrepolyder(x[:,4], N[4])    
        for i = 1:n
            p[i,:] = kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:])
            dp[i,1,:] = kron(kron(kron(dpx1[i,:], p2[i,:]), p3[i,:]), p4[i,:])
            dp[i,2,:] = kron(kron(kron(p1[i,:], dpx2[i,:]), p3[i,:]), p4[i,:])
            dp[i,3,:] = kron(kron(kron(p1[i,:], p2[i,:]), dpx3[i,:]), p4[i,:])
            dp[i,4,:] = kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), dpx4[i,:])
        end    
    elseif dim==5
        p3, dpx3 = legendrepolyder(x[:,3], N[3])    
        p4, dpx4 = legendrepolyder(x[:,4], N[4])   
        p5, dpx5 = legendrepolyder(x[:,5], N[5])    
        for i = 1:n
            p[i,:] = kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:])
            dp[i,1,:] = kron(kron(kron(kron(dpx1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:])
            dp[i,2,:] = kron(kron(kron(kron(p1[i,:], dpx2[i,:]), p3[i,:]), p4[i,:]), p5[i,:])
            dp[i,3,:] = kron(kron(kron(kron(p1[i,:], p2[i,:]), dpx3[i,:]), p4[i,:]), p5[i,:])
            dp[i,4,:] = kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), dpx4[i,:]), p5[i,:])
            dp[i,5,:] = kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), dpx5[i,:])
        end          
    elseif dim==6
        p3, dpx3 = legendrepolyder(x[:,3], N[3])    
        p4, dpx4 = legendrepolyder(x[:,4], N[4])   
        p5, dpx5 = legendrepolyder(x[:,5], N[5])    
        p6, dpx6 = legendrepolyder(x[:,6], N[6])    
        for i = 1:n
            p[i,:] = kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:])
            dp[i,1,:] = kron(kron(kron(kron(kron(dpx1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:])
            dp[i,2,:] = kron(kron(kron(kron(kron(p1[i,:], dpx2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:])
            dp[i,3,:] = kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), dpx3[i,:]), p4[i,:]), p5[i,:]), p6[i,:])
            dp[i,4,:] = kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), dpx4[i,:]), p5[i,:]), p6[i,:])
            dp[i,5,:] = kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), dpx5[i,:]), p6[i,:])
            dp[i,6,:] = kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), dpx6[i,:])
        end          
    elseif dim==7
        p3, dpx3 = legendrepolyder(x[:,3], N[3])    
        p4, dpx4 = legendrepolyder(x[:,4], N[4])   
        p5, dpx5 = legendrepolyder(x[:,5], N[5])    
        p6, dpx6 = legendrepolyder(x[:,6], N[6])    
        p6, dpx6 = legendrepolyder(x[:,7], N[7])    
        for i = 1:n
            p[i,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:]), p7[i,:])
            dp[i,1,:] = kron(kron(kron(kron(kron(kron(dpx1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:]), p7[i,:])
            dp[i,2,:] = kron(kron(kron(kron(kron(kron(p1[i,:], dpx2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:]), p7[i,:])
            dp[i,3,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), dpx3[i,:]), p4[i,:]), p5[i,:]), p6[i,:]), p7[i,:])
            dp[i,4,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), dpx4[i,:]), p5[i,:]), p6[i,:]), p7[i,:])
            dp[i,5,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), dpx5[i,:]), p6[i,:]), p7[i,:])
            dp[i,6,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), dpx6[i,:]), p7[i,:])
            dp[i,7,:] = kron(kron(kron(kron(kron(kron(p1[i,:], p2[i,:]), p3[i,:]), p4[i,:]), p5[i,:]), p6[i,:]), dpx7[i,:])
        end                  
    else
        error("Dimension must not exceed 7")  
    end

    return p, dp
end

function ref2dom(x::Vector{Float64}, ymin::Float64, ymax::Float64)

    y = ymin .+ (0.5*(ymax-ymin))*(x .+ 1) 
    return y
end

function ref2dom(x::Matrix{Float64}, ymin::Vector{Float64}, ymax::Vector{Float64})

    y = 0.0*x

    n = size(x,2)
    for i = 1:n
        y[:,i] = ymin[i] .+ (0.5*(ymax[i]-ymin[i]))*(x[:,i] .+ 1) 
    end

    return y
end

function tensornodes(ymin::Vector{Float64}, ymax::Vector{Float64}, N::Vector{Int64})

    x = tensornodes(N)
    if length(N) == 1
        y = ref2dom(x, ymin[1], ymax[1])    
        return reshape(y,(N[1]+1,1))
    end
    y = ref2dom(x, ymin, ymax)    

    return y
end

function dom2ref(y::Vector{Float64}, ymin::Float64, ymax::Float64)

    x = (2/(ymax-ymin))*(y .- ymin) .- 1

    return x
end

function dom2ref(y::Matrix{Float64}, ymin::Vector{Float64}, ymax::Vector{Float64})

    x = 0.0*y

    n = size(y,2)
    for i = 1:n
        x[:,i] = (2/(ymax[i] - ymin[i]))*(y[:,i] .- ymin[i]) .- 1
    end

    return x
end

