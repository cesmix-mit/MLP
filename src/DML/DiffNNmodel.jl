using Revise, Enzyme  
using Functors: @functor

mutable struct DiffNN{T}
    convweight::Array{T,3}       # weight matrices for DNN layers
    convbias::Array{T,2}         # bias vectors for DNN layers    
    outputweight::Array{T,1} # weight vector for output layer
    activationweight::Array{T,2} # weight vectors for activation functions 
    activationscale::Array{T,2}  # scale vectors for activation functions 
    f::Function          # activation function for DNN layers
    df::Function         # derivative of activation function
end

function DiffNN(num_inputs, num_layers, f, initW, initb, initc, initw, inits, T::DataType)
    # initialize weight matrices and vectors
    convweight = T.(initW(num_inputs, num_inputs, num_layers)) 
    convbias = T.(initb(num_inputs, num_layers)) 
    if typeof(initc) == Array{T,1}
        outputweight = initc
    else
        outputweight = T.(initc(num_inputs)) 
    end    
    activationweight = T.(initw(num_inputs, num_layers)) 
    activationscale = T.(inits(num_inputs, num_layers)) 
    df(x) = gradient(f, x)[1]  # derivative of activation function
    DiffNN(convweight, convbias, outputweight, activationweight, activationscale, f, df)
end

function DiffNN(num_inputs, num_layers, f, df, initW, initb, initc, initw, inits, T::DataType)
    # initialize weight matrices and vectors
    convweight = T.(initW(num_inputs, num_inputs, num_layers)) 
    convbias = T.(initb(num_inputs, num_layers)) 
    if typeof(initc) == Array{T,1}
        outputweight = initc
    else
        outputweight = T.(initW(num_inputs)) 
    end    
    activationweight = T.(initw(num_inputs, num_layers)) 
    activationscale = T.(inits(num_inputs, num_layers)) 
    DiffNN(convweight, convbias, outputweight, activationweight, activationscale, f, df)
end

function DiffNN(options::ModelOptions, f, df, T::DataType=Float64)    
    DiffNN(options.num_inputs, options.num_layers, f, df, options.initW, options.initb, options.initc, options.initw, options.inits, T)
end

function DiffNN(options::ModelOptions, initc, f, df, T::DataType=Float64)    
    DiffNN(options.num_inputs, options.num_layers, f, df, options.initW, options.initb, initc, options.initw, options.inits, T)
end

@functor DiffNN

function DiffNN_Output(convweight, convbias, outputweight, activationweight, activationscale, f, D)
        
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
    num_descriptors = length(D)    

    # Forward calculation of the energy 
    Dout = 1.0*D
    for i = 1:num_layers        
        c = convweight[:,:,i]*Dout + convbias[:,i]
        for j = 1:num_descriptors            
            Dout[j] = Dout[j] + f(c[j], activationweight[j,i], activationscale[j,i])
        end
    end
    e = sum(outputweight.*Dout)

    return e 
end

function DiffNN_Output(dml::DiffNN, D) 
    return DiffNN_Output(dml.convweight, dml.convbias, dml.outputweight, dml.activationweight, 
            dml.activationscale, dml.f, D) 
end
function (dml::DiffNN{T})(D::Array{T,1}) where {T<:Real}
    return DiffNN_Output(dml, D) 
end
function (dml::DiffNN{T})(D::Array{T,2}) where {T<:Real}
    num_configs = size(D,2)
    energy = zeros(num_configs)
    for i = 1:num_configs 
        energy[i] =  DiffNN_Output(dml, D[:,i])          
    end
    return energy
end

function DiffNN_ReverseDiff(convweight, convbias, outputweight, activationweight, activationscale, f, df, D, dDdR) 
                
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
        c = convweight[:,:,i]*Dout[:,i] + convbias[:,i]
        for j = 1:num_descriptors            
            Dout[j,i+1] = Dout[j,i] + f(c[j], activationweight[j,i], activationscale[j,i])
        end
    end
    e = sum(outputweight.*Dout[:,num_layers+1])

    # Reverse differentiation of the energy to calculate forces   
    dedD = 1.0*outputweight
    for i = num_layers:-1:1 
        c = convweight[:,:,i]*Dout[:,i] + convbias[:,i]
        for j = 1:num_descriptors            
            c[j] = df(c[j], activationweight[j,i], activationscale[j,i])
        end
        dedD = dedD + convweight[:,:,i]' * (c .* dedD)        
    end
    Forces = reshape(dedD,(1, num_descriptors))*dDdR

    return e, Forces, dedD 
end

function DiffNN_ReverseDiff(dml::DiffNN, D, dDdR) 
    return DiffNN_ReverseDiff(dml.convweight, dml.convbias, dml.outputweight, dml.activationweight, 
                dml.activationscale,dml.f, dml.df, D, dDdR) 
end
function (dml::DiffNN{T})(D::Array{T,1}, dDdR::Array{T,2}) where {T<:Real}
    return DiffNN_ReverseDiff(dml, D, dDdR) 
end
function (dml::DiffNN{T})(D::Array{T,2}, dDdR::Array{T,3}) where {T<:Real}
    num_configs = size(D,2)
    num_derivatives = size(dDdR,2)
    energy = zeros(num_configs)
    forces = zeros(num_derivatives, num_configs)
    for i = 1:num_configs 
        energy[i], forces[:,i], dedD =  DiffNN_ReverseDiff(dml, D[:,i], dDdR[:,:,i]) 
    end
    return energy, forces
end
function (dml::DiffNN{T})(x::Tuple) where {T<:Real}
    return dml(x...)
end

function DiffNN_LossFuncion(convweight, convbias, outputweight, activationweight, 
    activationscale, f, df, lossfunc, D, dDdR, E, F, beta) 
                
    num_descriptors, num_configs = size(D)    
    num_layers = size(convweight,3)
    num_derivatives = size(dDdR,2)
    F = reshape(F, (num_derivatives, num_configs))
    
    loss = 0.0 
    for n = 1:num_configs
        e, Forces,~ = DiffNN_ReverseDiff(convweight, convbias, outputweight, activationweight, 
                        activationscale, f, df, D[:,n], dDdR[:,:,n])                 
        loss = loss + lossfunc(e, Forces, E[n], F[:,n], beta, num_derivatives)
    end

    return loss
end

function DiffNN_LossFuncion(convweight, convbias, outputweight, activationweight, activationscale, f, lossfunc, D, E) 
                
    num_descriptors, num_configs = size(D)    
    num_layers = size(convweight,3)
    
    loss = 0.0 
    for n = 1:num_configs        
        e = DiffNN_Output(convweight, convbias, outputweight, activationweight, activationscale, f, D[:,n])        
        loss = loss + lossfunc(e, E[n])
    end

    return loss
end

function LossFunction(dml::DiffNN, options::ModelOptions, X, Y) where {T<:Real}
    loss = 0.0 
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        loss = DiffNN_LossFuncion(dml.convweight, dml.convbias, dml.outputweight, 
            dml.activationweight, dml.activationscale, dml.f, options.lossfunc, X[1], Y[1]) 
    else
        loss = DiffNN_LossFuncion(dml.convweight, dml.convbias, dml.outputweight, dml.activationweight, 
            dml.activationscale, dml.f, dml.df, options.lossfunc, X[1], X[2], Y[1], Y[2], options.β) 
    end
    return loss 
end

function DiffNN_LossFunction(convweight, convbias, outputweight, activationweight, activationscale, 
    f, df, lossfunc, D, dDdR, E, F, beta, num_descriptors, num_layers, num_derivatives, num_configs) 
                
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
                c[j] = f(c[j], activationweight[j,i], activationscale[j,i])
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
                c[j] = df(c[j], activationweight[j,i], activationscale[j,i])
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

function DiffNN_LossFunction(convweight, convbias, outputweight, activationweight, activationscale, f, lossfunc, D, E, 
    num_descriptors, num_layers, num_configs) 
                
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
                c[j] = f(c[j], activationweight[j,i], activationscale[j,i])
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

function DiffNN_LossGradient(dml::DiffNN, lossfunc, D::Array{T,2}, dDdR::Array{T,3}, E) where {T<:Real}
    num_descriptors, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_derivatives = size(dDdR,2)

    dW = zeros(num_descriptors, num_descriptors, num_layers)
    db = zeros(num_descriptors, num_layers)
    dc = zeros(num_descriptors)
    dw = zeros(num_descriptors, num_layers)
    ds = zeros(num_descriptors, num_layers)
    Enzyme.autodiff(DiffNN_LossFunction, Active, Duplicated(dml.convweight, dW), Duplicated(dml.convbias, db), 
    Duplicated(dml.outputweight, dc), Duplicated(dml.activationweight, dw), Duplicated(dml.activationscale, ds), 
    dml.f, lossfunc, D, E, num_descriptors, num_layers, num_configs, num_derivatives)

    return dW, db, dc, dw, ds
end

function DiffNN_LossGradient(dml::DiffNN, lossfunc, D::Array{T,2}, dDdR::Array{T,3}, E, F, beta) where {T<:Real}
    num_descriptors, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_derivatives = size(dDdR,2)
    F = reshape(F, (num_derivatives, num_configs))

    dW = zeros(num_descriptors, num_descriptors, num_layers)
    db = zeros(num_descriptors, num_layers)
    dc = zeros(num_descriptors)
    dw = zeros(num_descriptors, num_layers)
    ds = zeros(num_descriptors, num_layers)
    Enzyme.autodiff(DiffNN_LossFunction, Active, Duplicated(dml.convweight, dW), Duplicated(dml.convbias, db), 
    Duplicated(dml.outputweight, dc), Duplicated(dml.activationweight, dw), Duplicated(dml.activationscale, ds), 
    dml.f, dml.df, lossfunc, D, dDdR, E, F, beta, num_descriptors, num_layers, num_derivatives, num_configs)

    return dW, db, dc, dw, ds 
end

function LossGradient(dml::DiffNN, options::ModelOptions, X, Y) where {T<:Real}
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        dW, db, dc, dw, ds = DiffNN_LossGradient(dml, options.lossfunc, X[1], X[2], Y[1]) 
    else
        dW, db, dc, dw, ds = DiffNN_LossGradient(dml, options.lossfunc, X[1], X[2], Y[1], Y[2], options.β) 
    end
    grad = (dW, db, dc, dw, ds)
    return grad
end

function DiffNN_LossFunctionEnzyme(dml::DiffNN, lossfunc, D::Array{T,3}, E) where {T<:Real}
    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(dml.convweight,3)

    loss = DiffNN_LossFunction(dml.convweight, dml.convbias, dml.outputweight, dml.activationweight,
        dml.activationscale, dml.f, lossfunc, D, E, num_descriptors, num_layers, num_configs) 

    return loss         
end

function DiffNN_LossFunctionEnzyme(dml::DiffNN, lossfunc, D::Array{T,3}, dDdR::Array{T,4}, E, F, beta) where {T<:Real}

    num_descriptors, num_atoms, num_configs = size(D)    
    num_layers = size(dml.convweight,3)
    num_derivatives = size(dDdR,2)
    F = reshape(F, (num_derivatives, num_configs))
    
    loss = DiffNN_LossFunction(dml.convweight, dml.convbias, dml.outputweight, dml.activationweight,
            dml.activationscale,dml.f, dml.df, lossfunc, D, dDdR, E, F, beta, num_descriptors, 
            num_layers, num_derivatives, num_configs) 

    return loss
end

function LossFunctionEnzyme(dml::DiffNN, options::ModelOptions, X, Y) where {T<:Real}
    loss = 0.0 
    m = methods(options.lossfunc)[1]    
    if m.nargs<5
        loss = DiffNN_LossFunctionEnzyme(dml, options.lossfunc, X[1], Y[1]) 
    else
        loss = DiffNN_LossFunctionEnzyme(dml, options.lossfunc, X[1], X[2], Y[1], Y[2], options.β) 
    end
    return loss 
end


