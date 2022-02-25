using Preprocessing

mutable struct Data
    filename::String     
    fileformat::String 
    fileextension::String 
    fileindex::Array{Int32,1}
    num_atomtypes::Int32     
    num_descriptors::Int32 
    num_outputs::Int32     
    num_atoms::Array{Int32,1} 
    num_derivatives::Array{Int32,1} 
    num_configurations::Array{Int32,1}   
    randomize::Bool     
end

function savedata(Ae, Af, be, bf, nf, filename, batchsize, nfiles, npath)
    
    M = size(Ae,2)    
    nconfigs = length(nf)

    k = Array(batchsize:batchsize:nconfigs)
    num_forces = nf[1]
    split = 0
    if maximum(abs.(nf .- nf[1])) > 0        
        split = 1
        @warn "Configurations in the data set $npath have diffirent number of atoms"
        k = Array(1:1:nconfigs)
        num_forces = Int32.(1.0*nf)
    else
        if nconfigs <= batchsize
            k = [nconfigs]             
        elseif k[end] < nconfigs
            k = [k; nconfigs]
        end         
        num_forces = Int32.(ones(length(k))*nf[1])       
    end
    k = [0; k]    

    p = 0
    out = [];    
    num_forces_sum = cumsum([Int32(0); num_forces])
    for j = 1:(length(k)-1)
        nforces = num_forces[j]
        ei = (k[j]+1):(k[j+1])
        if split==0
            fi = (k[j]*nforces+1):(k[j+1]*nforces)
        else
            fi = (num_forces_sum[j]+1):(num_forces_sum[j+1])
            #display([ei[1] ei[end] fi[1] fi[end]])
        end

        J = k[j+1] - k[j]
        B = reshape(permutedims(reshape(Ae[ei,:], (J, M)), [2 1]), (M,J))
        D = reshape(permutedims(reshape(Af[fi,:], (nforces, J, M)), [3 1 2]), (M*nforces, J))

        E = reshape(be[ei], (1, J))        
        F = reshape(bf[fi], (nforces, J))
        data = [B; D; E; F];                    
        p = p + 1
        if p == 1
            out = [nforces J]
        else
            out = [out; [nforces J]]
        end
        numfile = nfiles + p 
        Preprocessing.savebin(filename * "$numfile.bin",data[:]);        
    end
    return out 
end

function DNNdata(datapath, descriptors, options, filename="trainbatch", batchsize=32,potential=nothing) 

    datapath = Preprocessing.deletenothing(datapath)
    nspecies = length(datapath[1].atomspecies)       
    fileformat = "bin"
    fileextension = "bin"
    randomize = false 
    num_outputs = 1        

    nfiles = 0 
    num_descriptors = 0
    num_atoms = []
    num_derivatives = []
    num_configurations = []
    for i = 1:length(datapath)
        # Ae : num_configurations by num_descriptors
        # Af : (dim * num_atoms * num_configurations) by num_descriptors
        Ae, Af, As, be, bf, bs, nf = globaldescriptors(datapath[i], descriptors, options, potential)     
        out = savedata(Ae, Af, be, bf, nf, filename, batchsize, nfiles, i)
        nfiles = nfiles + size(out,1)
        num_descriptors = size(Ae,2)            
        num_derivatives = [num_derivatives; Int32.(out[:,1])]    
        num_configurations = [num_configurations; Int32.(out[:,2])]      
        num_atoms = [num_atoms; Int32.(out[:,1]/3.0)]  
    end

    fileindex = Int32.(Array(1:nfiles))
    data = Data(filename, fileformat, fileextension, fileindex, nspecies, num_descriptors, 
            num_outputs, num_atoms, num_derivatives, num_configurations, randomize)
        
    return data 
end

function savedata(Ae, Af, be, bf, AA, PP, ne, nf, nspecies, nfiles, filename, batchsize, npath)
    p = 0    
    M = size(Ae,2)
    out = [];
    
    nconfigs = length(nf)

    natoms = ne[1]
    split = 0
    if maximum(abs.(ne .- natoms)) > 0
        @warn "Configurations in the data set $npath have diffirent number of atoms"
        split = 1
    end
    nAf = nf[1]
    if maximum(abs.(nf .- nAf)) > 0
        @warn "Configurations in the data set $npath have diffirent number of atoms"        
        split = 1
    end

    if split == 1        
        k = Array(1:1:nconfigs)
    else
        k = Array(batchsize:batchsize:nconfigs)
        if nconfigs <= batchsize
            k = [nconfigs]             
        elseif k[end] < nconfigs
            k = [k; nconfigs]
        end    
    end
    k = [0; k]

    nB_sum = cumsum([Int32(0); ne])
    nD_sum = cumsum([Int32(0); nf])
    nA_sum = cumsum([Int32(0); ne.*ne])
    nP_sum = cumsum([Int32(0); ne*nspecies])
    nF_sum = cumsum([Int32(0); Int32.(nf./ne)])
    for j = 1:(length(k)-1)               
        natoms = ne[j]
        nAf = nf[j]
        nforces = Int32(nAf/natoms)        
        if split==0            
            iB = (k[j]*natoms+1):(k[j+1]*natoms)
            iD = (k[j]*nAf+1):(k[j+1]*nAf)
            iA = (k[j]*natoms*natoms+1):(k[j+1]*natoms*natoms)    
            iP = (k[j]*natoms*nspecies+1):(k[j+1]*natoms*nspecies)    
            iE = (k[j]+1):(k[j+1])
            iF = (k[j]*nforces+1):(k[j+1]*nforces)    
        else
            iB = (nB_sum[j]+1):(nB_sum[j+1])
            iD = (nD_sum[j]+1):(nD_sum[j+1])             
            iA = (nA_sum[j]+1):(nA_sum[j+1])    
            iP = (nP_sum[j]+1):(nP_sum[j+1])    
            iE = (k[j]+1):(k[j+1])
            iF = (nF_sum[j]+1):(nF_sum[j+1])    
            #display([iB[1] iB[end] iD[1] iD[end]])           
        end

        J = k[j+1] - k[j] # number of configurations 
        B = reshape(permutedims(reshape(Ae[iB,:], (natoms, J, M)), [3 1 2]), (M*natoms,J))
        D = reshape(permutedims(reshape(Af[iD,:], (nforces, natoms, J, M)), [4 2 1 3]), (M*natoms*nforces, J))
        
        A = reshape(AA[iA], (natoms*natoms, J))                    
        P = reshape(PP[iP], (natoms*nspecies, J))            
        
        E = reshape(be[iE], (1, J))        
        F = reshape(bf[iF], (Int32(nforces), J))
        data = [B; D; A; P; E; F];                    
        p = p + 1
        if p == 1
            out = [natoms nforces J]
        else
            out = [out; [natoms nforces J]]
        end
        numfile = nfiles + p
        Preprocessing.savebin(filename * "$numfile.bin",data[:]);        
    end
    
    # nB = nB + M*natoms*nconfigs
    # nD = nD + M*natoms*nforces*nconfigs 
    # nE = nE + nconfigs
    # nA = nA + natoms*natoms*nconfigs 
    # nP = nP + natoms*nspecies*nconfigs 
    # nF = nF + nforces*nconfigs 
    
    return out 
end

function GNNdata(datapath, descriptors, options, filename="trainbatch", batchsize=32,potential=nothing) 

    datapath = Preprocessing.deletenothing(datapath)
    nspecies = length(datapath[1].atomspecies)       
    fileformat = "bin"
    fileextension = "bin"
    randomize = false 
    num_outputs = 1        

    nfiles = 0 
    num_descriptors = 0
    num_atoms = []
    num_derivatives = []
    num_configurations = []
    for i = 1:length(datapath)
        # Ae : (num_atoms * num_configurations) by num_descriptors
        # Af : (dim * num_atoms * num_atoms * num_configurations) by num_descriptors
        Ae, Af, be, bf, A, P, ne, nf = peratomdescriptors(datapath[i], descriptors, options, potential)     
        # display(size(Ae))
        # display(size(Af))
        out = savedata(Ae, Af, be, bf, A, P, ne, nf, nspecies, nfiles, filename, batchsize, i)
        nfiles = nfiles + size(out,1)
        num_descriptors = size(Ae,2)    
        num_atoms = [num_atoms; Int32.(out[:,1])]
        num_derivatives = [num_derivatives; Int32.(out[:,2])]    
        num_configurations = [num_configurations; Int32.(out[:,3])]        
    end
    
    fileindex = Int32.(Array(1:nfiles))
    data = Data(filename, fileformat, fileextension, fileindex, nspecies, num_descriptors, 
            num_outputs, num_atoms, num_derivatives, num_configurations, randomize)
        
    return data 
end



# # M : number of input global descriptors 
# # P : number of output global descriptors 
# # K : number of partial derivatives of the global descriptors 
# # J : batch size  
# # B : M by J          -> Input global descriptors 
# # D : M by K by J     -> Partial derivatives of input global descriptors 
# # E : P by J          -> Output global descriptors   
# # F : P by K by J     -> Partial derivatives of output global descriptors
# function getdata(data, M::Int32, P::Int32, K::Int32, T::DataType = Float64) 
#     nrB = M
#     nrD = M*K    
#     nrE = P
#     nrF = P*K

#     sz = size(data)
#     if length(sz)==1
#         I = Int32(nrB + nrD + nrE + nrF) 
#         J = Int32(sz[1]/I) 
#         data = reshape(data,(I,J))
#     elseif length(sz)==2
#         I, J = size(data)
#     else
#         error("data is not readable because its size is incorrect.")
#     end

#     if abs(I - (nrB+nrD+nrE+nrF)) > 0 
#         error("data is not readable because its size is incorrect.")
#     end

#     # X = (B, D)
#     X = (T.(reshape(data[1:nrB,:], (M, J))), 
#          T.(reshape(data[(nrB+1):(nrB+nrD),:], (M, K, J))), 
#         )

#     # Y = (E, F)         
#     Y = (T.(reshape(data[(nrB+nrD+1):(nrB+nrD+nrE),:], (P, J))), 
#          T.(reshape(data[(nrB+nrD+nrE+1):(nrB+nrD+nrE+nrF),:], (P, K, J)))
#          )

#     return X, Y  
# end

# # N : number of nodes/atoms 
# # M : number of peratom descriptors 
# # P : number of output global descriptors 
# # K : number of partial derivatives of the peratom descriptors 
# # Q : number of species 
# # J : batch size  
# # B : M by N by J          -> Input descriptor matrices    
# # D : K by M by N by J     -> Partial derivatives of input descriptor matrices 
# # A : N by N by J          -> Adjacency or Laplacian matrices  
# # P : N by Q by J          -> Atom type vectors 
# # E : P by J               -> Output global descriptor vectors     
# # F : P by K by J          -> Partial derivatives of output global descriptor vectors   
# function getdata(data, M::Int32, P::Int32, N::Int32, K::Int32, Q::Int32, T::DataType = Float64) 
#     nrB = M*N
#     nrD = M*K*N
#     nrA = N*N
#     nrP = N*Q
#     nrE = P
#     nrF = P*K

#     sz = size(data)
#     if length(sz)==1
#         I = Int32(nrB + nrD + nrA + nrP + nrE + nrF) 
#         J = Int32(sz[1]/I) 
#         data = reshape(data,(I,J))
#     elseif length(sz)==2
#         I, J = size(data)
#     else
#         error("data is not readable because its size is incorrect.")
#     end

#     if abs(I - (nrB+nrD+nrA+nrP+nrE+nrF)) > 0 
#         error("data is not readable because its size is incorrect.")
#     end

#     # X = (B, D, A, P)
#     X = (T.(reshape(data[1:nrB,:], (M, N, J))), 
#          T.(reshape(data[(nrB+1):(nrB+nrD),:], (K, M, N, J))), 
#          T.(reshape(data[(nrB+nrD+1):(nrB+nrD+nrA),:], (N, N, J))), 
#          T.(reshape(data[(nrB+nrD+nrA+1):(nrB+nrD+nrA+nrP),:], (N, Q, J))))

#     # Y = (E, F)         
#     Y = (T.(reshape(data[(nrB+nrD+nrA+nrP+1):(nrB+nrD+nrA+nrP+nrE),:], (P, J))), 
#          T.(reshape(data[(nrB+nrD+nrA+nrP+nrE+1):(nrB+nrD+nrA+nrP+nrE+nrF),:], (P, K, J))))

#     return X, Y  
# end

# function loaddata(data::Data, index, batchsize, T::DataType = Float64) 
         
#     if index > length(data.fileindex)
#         error("data is not available.")
#     end

#     filename = data.filename     
#     fileextension = data.fileextension
#     fileformat = lowercase(data.fileformat)    
#     fileindex = data.fileindex[index]
    
#     if fileformat == "text"
#         dataarray = Main.DelimitedFiles.readdlm(filename * string(fileindex));
#     else
#         dataarray = reinterpret(Float64,read(filename * string(fileindex) * "." * fileextension));
#     end

#     num_configs = data.num_configurations[index]
#     num_rows = Int32.(length(dataarray)/num_configs)
#     if (length(dataarray) % num_configs) !== 0
#         error("data is not readable because its size is incorrect.")
#     end 

#     if data.randomize == 1
#         ind = randperm(num_configs)
#     else
#         ind = Array(1:num_configs)
#     end
#     batch = min(batchsize, num_configs)    
#     dataarray = reshape(dataarray, (num_rows, num_configs))
#     dataarray = dataarray[:,ind[1:batch]]

#     M = data.num_descriptors
#     P = data.num_outputs
#     N = data.num_atoms[index]
#     K = data.num_derivatives[index]
#     Q = data.num_atomtypes
#     nr1 = M + M*K + P + P*K
#     nr2 = M*N + K*M*N + N*N + N*Q + P + P*K                    

#     if num_rows == nr1
#         return getdata(dataarray, M, P, K, T)  
#     elseif num_rows == nr2
#         return getdata(dataarray, M, P, N, K, Q, T)  
#     else
#         error("data is not readable because its size is incorrect.")        
#     end
# end

# function MLPoutputcoeff(data::Data, T::DataType = Float64) 
#     M = data.num_descriptors
#     A = zeros(M, M)
#     R = zeros(M)    
#     for i = 1:length(data.fileindex)
#         X, Y = loaddata(data, i, 100000000, T)     
#         AE = X[1]
#         J = data.num_configurations[i]
#         K = data.num_derivatives[i]        
#         AF = reshape(X[2], (M, K*J))
#         A = A + AE*(AE') + (1.0/K)*AF*(AF');
#         R = R + AE*(Y[1][:]) + (1.0/K)*AF*(Y[2][:]);
#     end
#     coeff = A\R         
#     return coeff
#     # for i = 1:length(data.fileindex)
#     #     X, Y = loaddata(data, i, 1000000, T)     
#     #     AE = X[1]
#     #     J = data.num_configurations[i]
#     #     K = data.num_derivatives[i]        
#     #     AF = reshape(X[2], (M, K*J))        
#     #     e = (AE')*coeff[:] - Y[1][:]
#     #     display(e)
#     #     e = (AF')*coeff[:] - Y[2][:]
#     #     display(e)
#     # end    
# end


