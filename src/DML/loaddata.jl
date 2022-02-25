using Preprocessing

# M : number of input global descriptors 
# P : number of output global descriptors 
# K : number of partial derivatives of the global descriptors 
# J : batch size  
# B : M by J          -> Input global descriptors 
# D : M by K by J     -> Partial derivatives of input global descriptors 
# E : P by J          -> Output global descriptors   
# F : P by K by J     -> Partial derivatives of output global descriptors
function getDNNdata(data, M::Int32, P::Int32, K::Int32, T::DataType = Float64) 
    nrB = M
    nrD = M*K    
    nrE = P
    nrF = P*K

    sz = size(data)
    if length(sz)==1
        I = Int32(nrB + nrD + nrE + nrF) 
        J = Int32(sz[1]/I) 
        data = reshape(data,(I,J))
    elseif length(sz)==2
        I, J = size(data)
    else
        error("data is not readable because its size is incorrect.")
    end

    if abs(I - (nrB+nrD+nrE+nrF)) > 0 
        error("data is not readable because its size is incorrect.")
    end

    # X = (B, D)
    X = (T.(reshape(data[1:nrB,:], (M, J))), 
         T.(reshape(data[(nrB+1):(nrB+nrD),:], (M, K, J))), 
        )

    # Y = (E, F)         
    Y = (T.(reshape(data[(nrB+nrD+1):(nrB+nrD+nrE),:], (P, J))), 
         T.(reshape(data[(nrB+nrD+nrE+1):(nrB+nrD+nrE+nrF),:], (P, K, J)))
         )

    return X, Y  
end

# N : number of nodes/atoms 
# M : number of input node/atom descriptors 
# P : number of output global descriptors 
# K : number of partial derivatives of the node/atom descriptors 
# J : batch size  
# B : M by N by J          -> Input descriptor matrices    
# D : M by K by N by J     -> Partial derivatives of input descriptor matrices 
# A : N by N by J          -> Adjacency or Laplacian matrices  
# Z : N by J               -> Atom type vectors 
# E : P by J               -> Output global descriptor vectors     
# F : P by K by J          -> Partial derivatives of output global descriptor vectors   
function getGNNdata(data, M::Int32, P::Int32, N::Int32, K::Int32, Q::Int32, T::DataType = Float64) 
    nrB = M*N
    nrD = M*K*N
    nrA = N*N
    nrP = N*Q
    nrE = P
    nrF = P*K

    sz = size(data)
    if length(sz)==1
        I = Int32(nrB + nrD + nrA + nrP + nrE + nrF) 
        J = Int32(sz[1]/I) 
        data = reshape(data,(I,J))
    elseif length(sz)==2
        I, J = size(data)
    else
        error("data is not readable because its size is incorrect.")
    end

    if abs(I - (nrB+nrD+nrA+nrP+nrE+nrF)) > 0 
        error("data is not readable because its size is incorrect.")
    end

    # X = (B, D, A, P)
    X = (T.(reshape(data[1:nrB,:], (M, N, J))), 
         T.(reshape(data[(nrB+1):(nrB+nrD),:], (M, N, K, J))), 
         T.(reshape(data[(nrB+nrD+1):(nrB+nrD+nrA),:], (N, N, J))), 
         T.(reshape(data[(nrB+nrD+nrA+1):(nrB+nrD+nrA+nrP),:], (N, Q, J))))

    # Y = (E, F)         
    Y = (T.(reshape(data[(nrB+nrD+nrA+nrP+1):(nrB+nrD+nrA+nrP+nrE),:], (P, J))), 
         T.(reshape(data[(nrB+nrD+nrA+nrP+nrE+1):(nrB+nrD+nrA+nrP+nrE+nrF),:], (P, K, J))))

    return X, Y  
end

function LoadData(data, index, batchsize, T::DataType = Float64) 
         
    if index > length(data.fileindex)
        error("data is not available.")
    end

    filename = data.filename     
    fileextension = data.fileextension
    fileformat = lowercase(data.fileformat)    
    fileindex = data.fileindex[index]
    
    if fileformat == "text"
        dataarray = Main.DelimitedFiles.readdlm(filename * string(fileindex));
    else
        dataarray = reinterpret(Float64,read(filename * string(fileindex) * "." * fileextension));
    end

    num_configs = data.num_configurations[index]
    num_rows = Int32.(length(dataarray)/num_configs)
    if (length(dataarray) % num_configs) !== 0
        error("data is not readable because its size is incorrect.")
    end 

    if data.randomize == 1
        ind = randperm(num_configs)
    else
        ind = Array(1:num_configs)
    end
    batch = min(batchsize, num_configs)    
    dataarray = reshape(dataarray, (num_rows, num_configs))
    dataarray = dataarray[:,ind[1:batch]]

    M = data.num_descriptors
    P = data.num_outputs
    N = data.num_atoms[index]
    K = data.num_derivatives[index]
    Q = data.num_atomtypes
    nr1 = M + M*K + P + P*K    
    nr2 = M*N + K*M*N + N*N + N*Q + P + P*K                    

    if num_rows == nr1
        return getDNNdata(dataarray, M, P, K, T)  
    elseif num_rows == nr2
        return getGNNdata(dataarray, M, P, N, K, Q, T)  
    else
        error("data is not readable because its size is incorrect.")        
    end
end

function DNNoutputweight(data, options::ModelOptions, T::DataType = Float64) 
    M = data.num_descriptors
    A = zeros(M, M)
    R = zeros(M)    
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, 100000000, T)     
        AE = X[1]
        J = data.num_configurations[i]
        K = data.num_derivatives[i]   
        N = data.num_atoms[i]     
        AF = reshape(X[2], (M, K*J))
        if options.normalizeenergy == true
            A = A + AE*(AE') + AF*(AF');
            R = R + AE*(Y[1][:]) + AF*(Y[2][:]);                
        else
            A = A + (1.0/(N*N))*AE*(AE')  + AF*(AF');
            R = R + (1.0/(N*N))*AE*(Y[1][:]) + AF*(Y[2][:]);        
        end
    end
    outputweight = A\R             
    return outputweight
end

function DNNscaleweight(data, T::DataType = Float64) 
    M = data.num_descriptors
    nconfigs = 0
    Dmean = zeros(M) 
    for i = 1:length(data.fileindex)
        X,~ = LoadData(data, i, 100000000, T)             
        nconfigs = nconfigs + size(X[1],2)
        tm = sum(X[1],dims=2) 
        Dmean = Dmean + tm[:]
    end
    Dmean = abs.(Dmean/nconfigs)     
    Dmean = 1.0./Dmean 
    return Dmean 
end

function DNNscaleweight(data, initweight, T::DataType = Float64) 
    M = data.num_descriptors
    nconfigs = 0
    Dmean = zeros(M) 
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, 100000000, T)       
        AE = X[1]      
        for m = 1:M
            AE[m,:] = AE[m,:]*initweight[m]; 
        end        
        nconfigs = nconfigs + size(X[1],2)
        tm = sum(AE, dims=2) 
        Dmean = Dmean + tm[:]
    end
    Dmean = abs.(Dmean/nconfigs) 
    Dmean = (1.0./Dmean).*(initweight) 
    return Dmean 
end

function DNNscalingdescriptors(data, options::ModelOptions, T::DataType = Float64) 
    outputweight = DNNoutputweight(data, options, T) 
    initweight = DNNscaleweight(data, outputweight, T)     
    #initweight = DNNscaleweight(data, T)     
    #initweight = DNNoutputweight(data, options, T) 
    M = data.num_descriptors
    A = zeros(M, M)
    R = zeros(M)    
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, 100000000, T)     
        AE = X[1]
        J = data.num_configurations[i]
        K = data.num_derivatives[i]   
        N = data.num_atoms[i]     
        AF = reshape(X[2], (M, K*J))
        for m = 1:M
            AE[m,:] = AE[m,:]*initweight[m]; 
            AF[m,:] = AF[m,:]*initweight[m]; 
        end        
        if options.normalizeenergy == true
            A = A + AE*(AE') + AF*(AF');
            R = R + AE*(Y[1][:]) + AF*(Y[2][:]);
        else
            A = A + (1.0/(N*N))*AE*(AE') + AF*(AF');
            R = R + (1.0/(N*N))*AE*(Y[1][:]) + AF*(Y[2][:]);
        end
        dataarray = [AE; reshape(AF, (M*K,J)); reshape(Y[1],(1,J)); reshape(Y[2], (K,J))]
        Preprocessing.savebin(data.filename * "$i.bin",dataarray[:]);        
    end
    for i = 1:M
         A[i,i] = A[i,i]*(1.0 + options.regparam)
    end
    # show(stdout, "text/plain", A)
    postweight = A\R        
    show(stdout, "text/plain", [outputweight initweight postweight])
    return postweight, initweight
end

function GNNoutputweight(data, options::ModelOptions, T::DataType = Float64) 
    M = data.num_descriptors
    Q = data.num_atomtypes
    A = zeros(M*Q, M*Q)
    R = zeros(M*Q)    
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, 100000000, T)             
        J = data.num_configurations[i]
        N = data.num_atoms[i]
        K = data.num_derivatives[i]        
        AE = reshape(X[1], (M, N, J))
        AF = permutedims(reshape(X[2], (M, N, K, J)), [1 3 2 4])
        PP = reshape(X[4], (N, Q, J))
        EP = zeros(M,Q,J)
        FP = zeros(M,K,Q,J)
        for j = 1:J
            EP[:,:,j] = AE[:,:,j]*PP[:,:,j]     
            FP[:,:,:,j] = reshape(reshape(AF[:,:,:,j], (M*K, N))*PP[:,:,j], (M,K,Q))            
        end
        EP = reshape(EP, (M*Q, J))
        FP = reshape(permutedims(FP, [1 2 3 4]), (M*Q, K*J))
        
        if options.normalizeenergy == true
            A = A + EP*(EP') + FP*(FP');
            R = R + EP*(Y[1][:]) + FP*(Y[2][:]);        
        else
            A = A + (1.0/(N*N))*EP*(EP') + FP*(FP');
            R = R + (1.0/(N*N))*EP*(Y[1][:]) + FP*(Y[2][:]);        
        end
        # if i == 24
        #     display([M N K J])
        #     savebin("GNN_A.bin", A[:])
        #     savebin("GNN_R.bin", R[:])        
        # end
    end
    outputweight = A\R         
    return outputweight
end

function GNNscalingdescriptors(data, options::ModelOptions, T::DataType = Float64) 
    initweight = GNNoutputweight(data, options, T) 
    M = data.num_descriptors
    P = data.num_outputs
    Q = data.num_atomtypes
    A = zeros(M*Q, M*Q)
    R = zeros(M*Q)    
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, 100000000, T)             
        J = data.num_configurations[i]
        K = data.num_derivatives[i]   
        N = data.num_atoms[i]     
        AE = reshape(X[1], (M, N, J))
        AF = reshape(X[2], (M, N, K, J))        
        for m = 1:M
            AE[m,:,:] = AE[m,:,:]*initweight[m]; 
            AF[m,:,:,:] = AF[m,:,:,:]*initweight[m]; 
        end        
        
        dataarray = [reshape(AE, (M*N,J)); reshape(AF, (M*N*K,J)); reshape(X[3],(N*N,J)); reshape(X[4],(N*Q,J)); reshape(Y[1],(P,J)); reshape(Y[2],(P*K,J))]
        Preprocessing.savebin(data.filename * "$i.bin",dataarray[:]);        

        AF = permutedims(AF, [1 3 2 4])
        PP = reshape(X[4], (N, Q, J))
        EP = zeros(M,Q,J)
        FP = zeros(M,K,Q,J)
        for j = 1:J
            EP[:,:,j] = AE[:,:,j]*PP[:,:,j]     
            FP[:,:,:,j] = reshape(reshape(AF[:,:,:,j], (M*K, N))*PP[:,:,j], (M,K,Q))            
        end
        EP = reshape(EP, (M*Q, J))
        FP = reshape(permutedims(FP, [1 2 3 4]), (M*Q, K*J))

        if options.normalizeenergy == true
            A = A + EP*(EP') + FP*(FP');
            R = R + EP*(Y[1][:]) + FP*(Y[2][:]);                
        else
            A = A + (1.0/(N*N))*EP*(EP') + FP*(FP');
            R = R + (1.0/(N*N))*EP*(Y[1][:]) + FP*(Y[2][:]);                        
        end
    end
    postweight = A\R         
    return postweight, initweight
end
