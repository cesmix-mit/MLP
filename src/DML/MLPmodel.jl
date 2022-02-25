using Revise, Enzyme  
using Flux: @functor

mutable struct BaseMLP{T}
    convweight::Array{T,3}   # weight matrices for MLP layers
    convbias::Array{T,2}     # bias for MLP layers
    outputweight::Array{T,1} # weight vector for output layer
    f::Function          # activation function for MLP layers
    df::Function         # derivative of activation function
end

function BaseMLP(num_inputs, num_layers, f, initW, initb, initc, T::DataType)
    # initialize weight matrices and vectors
    convweight = T.(initW(num_inputs, num_inputs, num_layers)) 
    convbias = T.(initb(num_inputs, num_layers)) 
    if typeof(initc) == Array{T,1}
        outputweight = initc
    else
        outputweight = T.(initW(num_inputs)) 
    end    
    df(x) = gradient(f, x)[1]  # derivative of activation function
    BaseMLP(convweight, convbias, outputweight, f, df)
end

function BaseMLP(num_inputs, num_layers, f, df, initW, initb, initc, T::DataType)
    # initialize weight matrices and vectors
    convweight = T.(initW(num_inputs, num_inputs, num_layers)) 
    convbias = T.(initb(num_inputs, num_layers)) 
    if typeof(initc) == Array{T,1}
        outputweight = initc
    else
        outputweight = T.(initW(num_inputs)) 
    end    
    BaseMLP(convweight, convbias, outputweight, f, df)
end

function BaseMLP(options::ModelOptions, f, df, T::DataType=Float64)    
    BaseMLP(options.num_inputs, options.num_layers, f, df, options.initW, options.initb, options.initc, T)
end

function BaseMLP(options::ModelOptions, initc, f, df, T::DataType=Float64)    
    BaseMLP(options.num_inputs, options.num_layers, f, df, options.initW, options.initb, initc, T)
end

@functor BaseMLP

function MLP_Output(convweight, convbias, outputweight, f, D)
        
    if (length(size(convweight))) != 3
        error("The convolution weight must be a three-dimensional array.")
    end
    if size(convweight,1) != size(convweight,2)
        error("The convolution weight must be a three-dimensional array of (num_descriptors, num_descriptors, num_layers).")
    end
    if length(D) != length(outputweight)
        error("the length of the output weight vector is invalid.")
    end

    num_layers = size(convweight,3)

    # Forward calculation of the energy 
    Dout = 1.0*D
    for i = 1:num_layers
        Dout = Dout + f.(convweight[:,:,i]*Dout + convbias[:,i])
    end
    e = sum(outputweight.*Dout)

    return e 
end

function MLP_Output(mlp::BaseMLP, D) 
    return MLP_Output(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, D) 
end
function (mlp::BaseMLP{T})(D::Array{T,1}) where {T<:Real}
    return MLP_Output(mlp, D) 
end
function (mlp::BaseMLP{T})(D::Array{T,2}) where {T<:Real}
    num_configs = size(D,2)
    energy = zeros(num_configs)
    for i = 1:num_configs 
        energy[i] =  MLP_Output(mlp, D[:,i])          
    end
    return energy
end

function MLP_ReverseDiff(convweight, convbias, outputweight, f, df, D, dDdR) 
                
    if (length(size(convweight))) != 3
        error("The convolution weight must be a three-dimensional array.")
    end
    if size(convweight,1) != size(convweight,2)
        error("The convolution weight must be a three-dimensional array of (num_descriptors, num_descriptors, num_layers).")
    end
    if length(D) != length(outputweight)
        error("the length of the output weight vector is invalid.")
    end

    num_descriptors = length(D)    
    num_layers = size(convweight,3)

    # Forward calculation of the energy 
    Dout = zeros(num_descriptors, num_layers+1)
    Dout[:,1] = D 
    for i = 1:num_layers
        Dout[:,i+1] = Dout[:,i] + f.(convweight[:,:,i]*Dout[:,i] + convbias[:,i])
        # if i==1
        #     Dout[:,i] = D + f.(convweight[:,:,i]*D + convbias[:,i])
        # else
        #     Dout[:,i] = Dout[:,i-1] + f.(convweight[:,:,i]*Dout[:,i-1] + convbias[:,i])
        # end
    end
    e = sum(outputweight.*Dout[:,num_layers+1])

    # Reverse differentiation of the energy to calculate forces   
    dedD = 1.0*outputweight
    for i = num_layers:-1:1 
        dedD = dedD + convweight[:,:,i]' * (df.(convweight[:,:,i]*Dout[:,i] + convbias[:,i]) .* dedD)
        # if i == 1
        #     #Q = repeat(dedD, outer = [1, num_descriptors])
        #     #dedD = dedD + (convweight[:,:,i] .* dedD)' * df.(convweight[:,:,i]*D)
        #     #Q = repeat(reshape(df.(convweight[:,:,i]*D),(1,num_descriptors)), outer = [num_descriptors,1])
        #     dedD = dedD + (df.(convweight[:,:,i]*D + convbias[:,i]) .* convweight[:,:,i])' * dedD 
        # else            
        #     #dedD = dedD + (convweight[:,:,i] .* dedD)' * df.(convweight[:,:,i]*Dout[:,i-1])
        #     dedD = dedD + (df.(convweight[:,:,i]*Dout[:,i-1] + convbias[:,i]) .* convweight[:,:,i])' * dedD
        # end
    end
    Forces = reshape(dedD,(1, num_descriptors))*dDdR

    return e, Forces, dedD 
end

function MLP_ReverseDiff(mlp::BaseMLP, D, dDdR) 
    return MLP_ReverseDiff(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, mlp.df, D, dDdR) 
end
function (mlp::BaseMLP{T})(D::Array{T,1}, dDdR::Array{T,2}) where {T<:Real}
    return MLP_ReverseDiff(mlp, D, dDdR) 
end
function (mlp::BaseMLP{T})(D::Array{T,2}, dDdR::Array{T,3}) where {T<:Real}
    num_configs = size(D,2)
    num_derivatives = size(dDdR,2)
    energy = zeros(num_configs)
    forces = zeros(num_derivatives, num_configs)
    for i = 1:num_configs 
        energy[i], forces[:,i], dedD =  MLP_ReverseDiff(mlp, D[:,i], dDdR[:,:,i]) 
    end
    return energy, forces
end
function (mlp::BaseMLP{T})(x::Tuple) where {T<:Real}
    return mlp(x...)
end

function MLP_LossFuncion(convweight, convbias, outputweight, f, df, lossfunc, D, dDdR, E, F, beta) 
                
    num_descriptors, num_configs = size(D)    
    num_layers = size(convweight,3)
    num_derivatives = size(dDdR,2)
    F = reshape(F, (num_derivatives, num_configs))
    #num_derivatives = size(F,1)
    
    loss = 0.0 
    for n = 1:num_configs
        # # Forward calculation of the energy 
        # Dout = zeros(num_descriptors, num_layers+1)
        # Dout[:,1] = 1.0*D[:,n] 
        # for i = 1:num_layers
        #     Dout[:,i+1] = Dout[:,i] + f.(convweight[:,:,i]*Dout[:,i] + convbias[:,i])
        # end
        # e = sum(outputweight.*Dout[:,num_layers+1])

        # # Reverse differentiation of the energy to calculate forces   
        # dedD = 1.0*outputweight
        # for i = num_layers:-1:1         
        #     dedD = dedD + convweight[:,:,i]' * (df.(convweight[:,:,i]*Dout[:,i] + convbias[:,i]) .* dedD)
        # end
        # Forces = reshape(dedD,(1, num_descriptors))*dDdR[:,:,n]
        
        e, Forces,~ = MLP_ReverseDiff(convweight, convbias, outputweight, f, df, D[:,n], dDdR[:,:,n])                 
        loss = loss + lossfunc(e, Forces, E[n], F[:,n], beta, num_derivatives)
    end

    return loss
end

function MLP_LossFuncion(convweight, convbias, outputweight, f, lossfunc, D, E) 
                
    num_descriptors, num_configs = size(D)    
    num_layers = size(convweight,3)
    
    loss = 0.0 
    for n = 1:num_configs        
        # Forward calculation of the energy 
        # Dout = zeros(num_descriptors, num_layers+1)
        # Dout[:,1] = 1.0*D[:,n] 
        # for i = 1:num_layers
        #     Dout[:,i+1] = Dout[:,i] + f.(convweight[:,:,i]*Dout[:,i] + convbias[:,i])
        # end
        # e = sum(outputweight.*Dout[:,num_layers+1])        
        e = MLP_Output(convweight, convbias, outputweight, f, D[:,n])        
        loss = loss + lossfunc(e, E[n])
    end

    return loss
end

function LossFunction(mlp::BaseMLP, lossfunc, X, Y) where {T<:Real}
    loss = MLP_LossFuncion(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, lossfunc, X[1], Y[1]) 
    return loss 
end

function LossFunction(mlp::BaseMLP, lossfunc, beta, X, Y) where {T<:Real}
    loss = MLP_LossFuncion(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, 
        mlp.df, lossfunc, X[1], X[2], Y[1], Y[2], beta) 
    return loss 
end

function LossFunction(mlp::BaseMLP, options::ModelOptions, X, Y) where {T<:Real}
    loss = 0.0 
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        loss = MLP_LossFuncion(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, options.lossfunc, X[1], Y[1]) 
    else
        loss = MLP_LossFuncion(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, 
            mlp.df, options.lossfunc, X[1], X[2], Y[1], Y[2], options.β) 
    end
    return loss 
end

function MLP_LossFunction(convweight, convbias, outputweight, f, df, lossfunc, D, dDdR, E, F, beta, 
    num_descriptors, num_layers, num_derivatives, num_configs) 
                
    # Forward calculation of the energy 
    c = zeros(num_descriptors)
    Dout = zeros(num_descriptors, num_layers+1)
    dedD = zeros(num_descriptors)    
    dC = zeros(num_descriptors)    

    Energy = zeros(num_configs)
    Forces = zeros(num_derivatives,num_configs)

    loss = 0.0
    for n = 1:num_configs
        # Forward calculation of the energy 
        for i = 1:num_descriptors
            Dout[i,1] = D[i,n]
        end
        for i = 1:num_layers                
            for j = 1:num_descriptors
                c[j] = convbias[j,i]
                for k = 1:num_descriptors
                    c[j] = c[j] + convweight[j,k,i]*Dout[k,i]
                end            
                c[j] = f(c[j])
            end
            Dout[:,i+1] = Dout[:,i] + c
        end
        for i = 1:num_descriptors
            Energy[n] = Energy[n] + outputweight[i]*Dout[i,num_layers+1]
        end    

        # Reverse differentiation of the energy to calculate forces   
        for i = 1:num_descriptors
            dedD[i] = outputweight[i]
        end    
        for i = num_layers:-1:1 
            for j = 1:num_descriptors
                c[j] = convbias[j,i]
                for k = 1:num_descriptors
                    c[j] = c[j] + convweight[j,k,i]*Dout[k,i]
                end            
                c[j] = df(c[j])
            end
            for j = 1:num_descriptors           
                dC[j] = 0.0
                for k = 1:num_descriptors
                    dC[j] = dC[j] + c[k] * convweight[k,j,i] * dedD[k]
                end
            end
            dedD = dedD + dC
        end    

        for i = 1:num_derivatives
            for j = 1:num_descriptors
                Forces[i,n] = Forces[i,n] + dedD[j]*dDdR[j,i,n]
            end
        end

        loss = loss + lossfunc(Energy[n], Forces[:,n], E[n], F[:,n], beta, num_derivatives)
    end

    return loss
end

function MLP_LossFunction(convweight, convbias, outputweight, f, lossfunc, D, E, num_descriptors, num_layers, num_configs) 
                
    # Forward calculation of the energy 
    c = zeros(num_descriptors)
    Dout = zeros(num_descriptors, num_layers+1)
    Energy = zeros(num_configs)

    loss = 0.0
    for n = 1:num_configs
        # Forward calculation of the energy 
        for i = 1:num_descriptors
            Dout[i,1] = D[i,n]
        end
        for i = 1:num_layers                
            for j = 1:num_descriptors
                c[j] = convbias[j,i]
                for k = 1:num_descriptors
                    c[j] = c[j] + convweight[j,k,i]*Dout[k,i]
                end            
                c[j] = f(c[j])
            end
            Dout[:,i+1] = Dout[:,i] + c
        end
        for i = 1:num_descriptors
            Energy[n] = Energy[n] + outputweight[i]*Dout[i,num_layers+1]
        end    
        loss = loss + lossfunc(Energy[n], E[n])
    end

    return loss
end

function MLP_LossGradient(mlp::BaseMLP, lossfunc, D::Array{T,2}, E) where {T<:Real}
    num_descriptors, num_configs = size(D)    
    num_layers = size(mlp.convweight,3)

    dW = zeros(num_descriptors, num_descriptors, num_layers)
    db = zeros(num_descriptors, num_layers)
    dc = zeros(num_descriptors)
    Enzyme.autodiff(MLP_LossFunction, Active, Duplicated(mlp.convweight, dW), Duplicated(mlp.convbias, db), 
    Duplicated(mlp.outputweight, dc), mlp.f, lossfunc, D, E, num_descriptors, num_layers, num_configs)

    return dW, db, dc
end

function MLP_LossGradient(mlp::BaseMLP, lossfunc, D::Array{T,2}, dDdR::Array{T,3}, E, F, beta) where {T<:Real}
    num_descriptors, num_configs = size(D)    
    num_layers = size(mlp.convweight,3)
    num_derivatives = size(dDdR,2)
    F = reshape(F, (num_derivatives, num_configs))

    dW = zeros(num_descriptors, num_descriptors, num_layers)
    db = zeros(num_descriptors, num_layers)
    dc = zeros(num_descriptors)
    Enzyme.autodiff(MLP_LossFunction, Active, Duplicated(mlp.convweight, dW), Duplicated(mlp.convbias, db), 
    Duplicated(mlp.outputweight, dc), mlp.f, mlp.df, lossfunc, D, dDdR, E, F, beta, num_descriptors, 
    num_layers, num_derivatives, num_configs)

    return dW, db, dc
end

function LossGradient(mlp::BaseMLP, lossfunc, X, Y) where {T<:Real}
    dW, db, dc = MLP_LossGradient(mlp, lossfunc, X[1], Y[1]) 
    grad = (dW, db, dc)
    return grad
end

function LossGradient(mlp::BaseMLP, lossfunc, beta, X, Y) where {T<:Real}
    dW, db, dc = MLP_LossGradient(mlp, lossfunc, X[1], X[2], Y[1], Y[2], beta) 
    grad = (dW, db, dc)
    return grad
end

function LossGradient(mlp::BaseMLP, options::ModelOptions, X, Y) where {T<:Real}
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        dW, db, dc = MLP_LossGradient(mlp, options.lossfunc, X[1], Y[1]) 
    else
        dW, db, dc = MLP_LossGradient(mlp, options.lossfunc, X[1], X[2], Y[1], Y[2], options.β) 
    end
    grad = (dW, db, dc)
    return grad
end

function MLP_LossFunctionEnzyme(mlp::BaseMLP, lossfunc, D::Array{T,3}, E) where {T<:Real}
    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(mlp.convweight,3)

    loss = MLP_LossFunction(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, lossfunc, D, E, 
            num_descriptors, num_layers, num_configs) 

    return loss         
end

function MLP_LossFunctionEnzyme(mlp::BaseMLP, lossfunc, D::Array{T,3}, dDdR::Array{T,4}, E, F, beta) where {T<:Real}

    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(mlp.convweight,3)
    num_derivatives = size(dDdR,2)
    F = reshape(F, (num_derivatives, num_configs))
    
    loss = MLP_LossFunction(mlp.convweight, mlp.convbias, mlp.outputweight, mlp.f, mlp.df, lossfunc, D, dDdR, E, F, beta, 
        num_descriptors, num_layers, num_derivatives, num_configs) 

    return loss
end

function LossFunctionEnzyme(mlp::BaseMLP, options::ModelOptions, X, Y) where {T<:Real}
    loss = 0.0 
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        loss = MLP_LossFunctionEnzyme(mlp, options.lossfunc, X[1], Y[1]) 
    else
        loss = MLP_LossFunctionEnzyme(mlp, options.lossfunc, X[1], X[2], Y[1], Y[2], options.β) 
    end
    return loss 
end


