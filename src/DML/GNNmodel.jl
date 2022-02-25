using Revise, Enzyme  
using Functors: @functor

mutable struct BaseGNN{T}
    convweight::Array{T,3}   # weight matrices for GNN layers
    convbias::Array{T,2}     # bias for GNN layers
    outputweight::Array{T,1} # weight vector for output layer
    f::Function          # activation function for GNN layers
    df::Function         # derivative of activation function
end

function BaseGNN(num_inputs, num_layers, f, initW, initb, initc, T::DataType)
    # initialize weight matrices and vectors
    convweight = T.(initW(num_inputs, num_inputs, num_layers)) 
    convbias = T.(initb(num_inputs, num_layers)) 
    if typeof(initc) == Array{T,1}
        outputweight = initc
    else
        outputweight = T.(initW(num_inputs)) 
    end    
    df(x) = gradient(f, x)[1]  # derivative of activation function
    BaseGNN(convweight, convbias, outputweight, f, df)
end

function BaseGNN(num_inputs, num_layers, f, df, initW, initb, initc, T::DataType)
    # initialize weight matrices and vectors
    convweight = T.(initW(num_inputs, num_inputs, num_layers)) 
    convbias = T.(initb(num_inputs, num_layers)) 
    if typeof(initc) == Array{T,1}
        outputweight = initc
    else
        outputweight = T.(initW(num_inputs)) 
    end    
    BaseGNN(convweight, convbias, outputweight, f, df)
end

function BaseGNN(options::ModelOptions, f, df, T::DataType=Float64)    
    BaseGNN(options.num_inputs, options.num_layers, f, df, options.initW, options.initb, options.initc, T)
end

function BaseGNN(options::ModelOptions, initc, f, df, T::DataType=Float64)    
    BaseGNN(options.num_inputs, options.num_layers, f, df, options.initW, options.initb, initc, T)
end

@functor BaseGNN

function GNN_Output(convweight, convbias, outputweight, f, D, A, P)
        
    if (length(size(convweight))) != 3
        error("The convolution weight must be a three-dimensional array.")
    end
    if size(convweight,1) != size(convweight,2)
        error("The convolution weight must be a three-dimensional array of (num_descriptors, num_descriptors, num_layers).")
    end

    num_layers = size(convweight,3)
    num_atoms = size(D,2)

    # Forward convolution layers to compute GNN descriptors
    Dout = 1.0*D
    for i = 1:num_layers                
        U = densesparsematmul(convweight[:,:,i]*Dout, A); 
        Dout = Dout + f.(U + repeat(convbias[:,i], outer=[1 num_atoms]))        
    end
    # Pooling GNN descriptors
    Cout = Dout*P
    # Energy 
    e = sum(outputweight.*Cout[:])

    return e 
end

function GNN_Output(dml::BaseGNN, D, A, P) 
    return GNN_Output(dml.convweight, dml.convbias, dml.outputweight, dml.f, D, A, P) 
end
function (dml::BaseGNN{T})(D::Array{T,2}, A::Array{T,2}, P::Array{T,2}) where {T<:Real}
    return GNN_Output(dml, D, A, P) 
end
function (dml::BaseGNN{T})(D::Array{T,3}, A::Array{T,3}, P::Array{T,3}) where {T<:Real}
    num_configs = size(D,3)
    energy = zeros(num_configs)
    for i = 1:num_configs 
        energy[i] =  GNN_Output(dml, D[:,:,i], A[:,:,i], P[:,:,i])          
    end
    return energy
end

function GNN_ReverseDiff(convweight, convbias, outputweight, f, df, D, A, P, dDdR) 
                
    if (length(size(convweight))) != 3
        error("The convolution weight must be a three-dimensional array.")
    end
    if size(convweight,1) != size(convweight,2)
        error("The convolution weight must be a three-dimensional array of (num_descriptors, num_descriptors, num_layers).")
    end

    num_descriptors, num_atoms = size(D)    
    num_layers = size(convweight,3)
    num_atomtypes = size(P,2) 

    # Forward calculation of the energy 
    Dout = zeros(num_descriptors, num_atoms, num_layers+1)
    Dout[:,:,1] = 1.0*D 
    for i = 1:num_layers
        U = densesparsematmul(convweight[:,:,i]*Dout[:,:,i], A); 
        Dout[:,:,i+1] = Dout[:,:,i] + f.(U + repeat(convbias[:,i], outer=[1 num_atoms]))
    end
    Cout = Dout[:,:,num_layers+1]*P
    e = sum(outputweight.*Cout[:])
    #e = sum(outputweight.*Dout[:,num_layers+1])

    # Reverse differentiation of the energy to calculate forces   
    dedD = reshape(outputweight,(num_descriptors,num_atomtypes)) * (P'); 
    for i = num_layers:-1:1 
        U = densesparsematmul(convweight[:,:,i]*Dout[:,:,i], A); 
        V = df.(U + repeat(convbias[:,i], outer=[1 num_atoms])).*dedD;    
        dedD = dedD + densesparsematmul(V, convweight[:,:,i], A);
        #U = convweight[:,:,i]*Dout[:,:,i]*A + repeat(convbias[:,i], outer=[1 num_atoms])        
        #dedD = dedD + df.(U) .* (convweight[:,:,i]' * dedD * A')        
    end
    num_derivatives = size(dDdR,3)
    Forces = reshape(dedD,(1, num_descriptors*num_atoms))*reshape(dDdR,(num_descriptors*num_atoms, num_derivatives))
    Forces = Forces[:]
    
    return e, Forces, dedD 
end

function GNN_ReverseDiff(dml::BaseGNN, D, A, P, dDdR) 
    return GNN_ReverseDiff(dml.convweight, dml.convbias, dml.outputweight, dml.f, dml.df, D, A, P, dDdR) 
end
function (dml::BaseGNN{T})(D::Array{T,2}, A::Array{T,2}, P::Array{T,2}, dDdR::Array{T,3}) where {T<:Real}
    return GNN_ReverseDiff(dml, D, A, P, dDdR) 
end
function (dml::BaseGNN{T})(D::Array{T,3}, A::Array{T,3}, P::Array{T,3}, dDdR::Array{T,4}) where {T<:Real}
    num_configs = size(D,3)
    num_derivatives = size(dDdR,3)
    energy = zeros(num_configs)
    forces = zeros(num_derivatives, num_configs)
    for i = 1:num_configs 
        energy[i], forces[:,i], dedD =  GNN_ReverseDiff(dml, D[:,:,i], A[:,:,i], P[:,:,i], dDdR[:,:,:,i]) 
    end
    return energy, forces
end
function (dml::BaseGNN{T})(x::Tuple) where {T<:Real}
    return dml(x...)
end

function GNN_LossFuncion(convweight, convbias, outputweight, f, lossfunc, D, A, P, E)                 
    num_configs = size(D,3)            
    loss = 0.0 
    for n = 1:num_configs
        e = GNN_Output(convweight, convbias, outputweight, f, D[:,:,n], A[:,:,n], P[:,:,n])        
        loss = loss + lossfunc(e, E[n])
    end
    return loss
end

function GNN_LossFuncion(convweight, convbias, outputweight, f, df, lossfunc, D, A, P, dDdR, E, F, beta)                 
    num_configs = size(D,3)           
    num_derivatives = size(dDdR,3)    
    F = reshape(F, (num_derivatives, num_configs))
    loss = 0.0 
    for n = 1:num_configs        
        e, Forces,~ = GNN_ReverseDiff(convweight, convbias, outputweight, f, df, D[:,:,n], A[:,:,n], P[:,:,n], dDdR[:,:,:,n])         
        loss = loss + lossfunc(e, Forces, E[n], F[:,n], beta, num_derivatives)
    end
    return loss
end

function LossFunction(dml::BaseGNN, lossfunc, X, Y) where {T<:Real}
    loss = GNN_LossFuncion(dml.convweight, dml.convbias, dml.outputweight, dml.f, 
        lossfunc, X[1], X[3], X[4], Y[1]) 
    return loss 
end
function LossFunction(dml::BaseGNN, lossfunc, beta, X, Y) where {T<:Real}
    loss = GNN_LossFuncion(dml.convweight, dml.convbias, dml.outputweight, dml.f, 
        dml.df, lossfunc, X[1], X[3], X[4], X[2], Y[1], Y[2], beta) 
    return loss 
end

function LossFunction(dml::BaseGNN, options::ModelOptions, X, Y) where {T<:Real}
    loss = 0.0 
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        loss = GNN_LossFuncion(dml.convweight, dml.convbias, dml.outputweight, dml.f, 
            options.lossfunc, X[1], X[3], X[4], Y[1]) 
    else
        loss = GNN_LossFuncion(dml.convweight, dml.convbias, dml.outputweight, dml.f, 
            dml.df, options.lossfunc, X[1], X[3], X[4], X[2], Y[1], Y[2], options.β) 
    end
    return loss 
end

function GNN_LossFunction(convweight, convbias, outputweight, f, lossfunc, D, A, P, E, 
    num_descriptors, num_layers, num_atoms, num_species, num_configs) 
                
    # Forward calculation of the energy 
    Bout = zeros(num_descriptors, num_atoms)
    Cout = zeros(num_descriptors, num_species)
    Dout = zeros(num_descriptors, num_atoms, num_layers+1)    
    U = zeros(num_descriptors, num_atoms)
    Energy = zeros(num_configs)

    loss = 0.0
    for c = 1:num_configs
        # Forward calculation of the energy 
        for i = 1:num_descriptors
            for j = 1:num_atoms 
                Dout[i,j,1] = D[i,j,c]
            end
        end
        for l = 1:num_layers               
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    Bout[i,j] = 0.0
                    for k = 1:num_descriptors
                        Bout[i,j] = Bout[i,j] + convweight[i,k,l]*Dout[k,j,l]
                    end
                end
            end
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    U[i,j] = 0.0
                end
            end
            for j = 1:num_atoms
                for k = 1:num_atoms
                    if abs(A[k,j,c]) > 0                        
                        for i = 1:num_descriptors   
                            U[i,j] = U[i,j] + Bout[i,k]*A[k,j,c];
                        end
                    end
                end
            end
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    Dout[i,j,l+1] = Dout[i,j,l] + f(U[i,j] + convbias[i,l])                             
                end
            end        
        end
        for i = 1:num_descriptors
            for j = 1:num_species 
                Cout[i,j] = 0.0
                for k = 1:num_atoms
                    Cout[i,j] = Cout[i,j] + Dout[i,k,num_layers+1]*P[k,j,c]
                end
            end
        end
        for i = 1:num_descriptors 
            for j = 1:num_species
                Energy[c] = Energy[c] + outputweight[i + num_descriptors*(j-1)]*Cout[i,j]
            end
        end    
        loss = loss + lossfunc(Energy[c], E[c])
    end

    return loss
end

function GNN_LossFunction(convweight, convbias, outputweight, f, df, lossfunc, D, A, P, dDdR, E, F, beta, 
    num_descriptors, num_layers, num_atoms, num_species, num_derivatives, num_configs) 
                
    # Forward calculation of the energy         
    Bout = zeros(num_descriptors, num_atoms)
    Cout = zeros(num_descriptors, num_species)
    Dout = zeros(num_descriptors, num_atoms, num_layers+1)    
    U = zeros(num_descriptors, num_atoms)
    V = zeros(num_descriptors, num_atoms)    
    W = zeros(num_descriptors, num_atoms)    
    dedD = zeros(num_descriptors, num_atoms)    

    Energy = zeros(num_configs)
    Forces = zeros(num_derivatives, num_configs)
    
    loss = 0.0
    for c = 1:num_configs        
        for i = 1:num_descriptors
            for j = 1:num_atoms 
                Dout[i,j,1] = D[i,j,c]
            end
        end
        for l = 1:num_layers               
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    Bout[i,j] = 0.0
                    for k = 1:num_descriptors
                        Bout[i,j] = Bout[i,j] + convweight[i,k,l]*Dout[k,j,l]
                    end
                end
            end
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    U[i,j] = 0.0
                end
            end
            for j = 1:num_atoms
                for k = 1:num_atoms
                    if abs(A[k,j,c]) > 0                        
                        for i = 1:num_descriptors   
                            U[i,j] = U[i,j] + Bout[i,k]*A[k,j,c];
                        end
                    end
                end
            end
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    Dout[i,j,l+1] = Dout[i,j,l] + f(U[i,j] + convbias[i,l])                             
                end
            end        
        end
        for i = 1:num_descriptors
            for j = 1:num_species 
                Cout[i,j] = 0.0
                for k = 1:num_atoms
                    Cout[i,j] = Cout[i,j] + Dout[i,k,num_layers+1]*P[k,j,c]
                end
            end
        end
        for i = 1:num_descriptors 
            for j = 1:num_species
                Energy[c] = Energy[c] + outputweight[i + num_descriptors*(j-1)]*Cout[i,j]
            end
        end    

        for i = 1:num_descriptors             
            for j = 1:num_atoms
                dedD[i,j] = 0.0
                for k = 1:num_species
                    dedD[i,j] = dedD[i,j] + outputweight[i + num_descriptors*(k-1)]*P[j,k,c]
                end
            end
        end    
        for l = num_layers:-1:1         
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    Bout[i,j] = 0.0
                    for k = 1:num_descriptors
                        Bout[i,j] = Bout[i,j] + convweight[i,k,l]*Dout[k,j,l]
                    end
                end
            end
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    U[i,j] = 0.0
                end
            end
            for j = 1:num_atoms
                for k = 1:num_atoms
                    if abs(A[k,j,c]) > 0                        
                        for i = 1:num_descriptors   
                            U[i,j] = U[i,j] + Bout[i,k]*A[k,j,c];
                        end
                    end
                end
            end
            for i = 1:num_descriptors
                for j = 1:num_atoms 
                    V[i,j] = dedD[i,j]*df(U[i,j] + convbias[i,l])                             
                end
            end        
            for m = 1:num_descriptors        
                for n = 1:num_atoms   
                    W[m,n] = 0.0
                end
            end
            for n = 1:num_atoms   
                for j = 1:num_atoms    
                    if abs(A[n,j,c]) > 0
                        for m = 1:num_descriptors        
                            for i = 1:num_descriptors
                                W[m,n] = W[m,n] + V[i,j]*convweight[i,m,l]*A[n,j,c];
                            end
                        end
                    end
                end
            end        
            for m = 1:num_descriptors        
                for n = 1:num_atoms   
                    dedD[m,n] = dedD[m,n] + W[m,n]
                end
            end
        end
        for i = 1:num_derivatives
            for k = 1:num_atoms                
                for j = 1:num_descriptors                
                    if abs(dDdR[j,k,i,c]) > 0
                        Forces[i,c] = Forces[i,c] + dedD[j,k]*dDdR[j,k,i,c]
                    end
                end
            end
        end
        loss = loss + lossfunc(Energy[c], Forces[:,c], E[c], F[:,c], beta, num_derivatives)
    end

    return loss
end

function GNN_LossGradient(dml::BaseGNN, lossfunc, D::Array{T,3}, A::Array{T,3}, P::Array{T,3}, E) where {T<:Real}
    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_species = size(P,2)

    dW = zeros(num_descriptors, num_descriptors, num_layers)
    db = zeros(num_descriptors, num_layers)
    dc = zeros(num_descriptors*num_species)
    Enzyme.autodiff(GNN_LossFunction, Active, Duplicated(dml.convweight, dW), Duplicated(dml.convbias, db), 
    Duplicated(dml.outputweight, dc), dml.f, lossfunc, D, A, P, E, num_descriptors, 
    num_layers, num_atoms, num_species, num_configs)

    return dW, db, dc
end

function GNN_LossGradient(dml::BaseGNN, lossfunc, D::Array{T,3}, A::Array{T,3}, P::Array{T,3}, 
    dDdR::Array{T,4}, E, F, beta) where {T<:Real}

    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_species = size(P,2)    
    num_derivatives = size(dDdR,3)
    F = reshape(F, (num_derivatives, num_configs))

    dW = zeros(num_descriptors, num_descriptors, num_layers)
    db = zeros(num_descriptors, num_layers)
    dc = zeros(num_descriptors*num_species)
    Enzyme.autodiff(GNN_LossFunction, Active, Duplicated(dml.convweight, dW), Duplicated(dml.convbias, db), 
    Duplicated(dml.outputweight, dc), dml.f, dml.df, lossfunc, D, A, P, dDdR, E, F, beta, num_descriptors, 
    num_layers, num_atoms, num_species, num_derivatives, num_configs)

    return dW, db, dc
end

function LossGradient(dml::BaseGNN, lossfunc, X, Y) where {T<:Real}
    dW, db, dc = GNN_LossGradient(dml, lossfunc, X[1], X[3], X[4], Y[1]) 
    grad = (dW, db, dc)
    return grad
end

function LossGradient(dml::BaseGNN, lossfunc, beta, X, Y) where {T<:Real}
    dW, db, dc = GNN_LossGradient(dml, lossfunc, X[1], X[3], X[4], X[2], Y[1], Y[2], beta) 
    grad = (dW, db, dc)
    return grad
end

function LossGradient(dml::BaseGNN, options::ModelOptions, X, Y) where {T<:Real}
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        dW, db, dc = GNN_LossGradient(dml, options.lossfunc, X[1], X[3], X[4], Y[1]) 
    else
        dW, db, dc = GNN_LossGradient(dml, options.lossfunc, X[1], X[3], X[4], X[2], Y[1], Y[2], options.β) 
    end
    grad = (dW, db, dc)
    return grad
end

function GNN_LossFunctionEnzyme(dml::BaseGNN, lossfunc, D::Array{T,3}, A::Array{T,3}, P::Array{T,3}, E) where {T<:Real}
    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_species = size(P,2)

    loss = GNN_LossFunction(dml.convweight, dml.convbias, dml.outputweight, dml.f, lossfunc, D, A, P, E, 
        num_descriptors, num_layers, num_atoms, num_species, num_configs) 

    return loss         
end

function GNN_LossFunctionEnzyme(dml::BaseGNN, lossfunc, D::Array{T,3}, A::Array{T,3}, P::Array{T,3}, 
    dDdR::Array{T,4}, E, F, beta) where {T<:Real}

    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_species = size(P,2)    
    num_derivatives = size(dDdR,3)
    F = reshape(F, (num_derivatives, num_configs))

    loss = GNN_LossFunction(dml.convweight, dml.convbias, dml.outputweight, dml.f, dml.df, lossfunc, D, A, P, dDdR, E, F, beta, 
        num_descriptors, num_layers, num_atoms, num_species, num_derivatives, num_configs) 

    return loss
end

function LossFunctionEnzyme(dml::BaseGNN, options::ModelOptions, X, Y) where {T<:Real}
    loss = 0.0 
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        loss = GNN_LossFunctionEnzyme(dml, options.lossfunc, X[1], X[3], X[4], Y[1]) 
    else
        loss = GNN_LossFunctionEnzyme(dml, options.lossfunc, X[1], X[3], X[4], X[2], Y[1], Y[2], options.β) 
    end
    return loss 
end

