function BuildModel(data, options::ModelOptions)            

    options.num_inputs = data.num_descriptors

    # Initialize DNN model 
    if options.modeltype == "dnn"        
        #outputweight = DNNoutputweight(data) 
        outputweight, initweight = DNNscalingdescriptors(data, options) 
        model = BaseDNN(options, outputweight, options.f, options.df) 
    elseif options.modeltype == "diffnn"        
        #outputweight = DNNoutputweight(data) 
        outputweight, initweight = DNNscalingdescriptors(data, options) 
        model = DiffNN(options, outputweight, options.f, options.df)     
    elseif options.modeltype == "gnn"
        #outputweight = GNNoutputweight(data)         
        outputweight, initweight = GNNscalingdescriptors(data, options) 
        model = BaseGNN(options, outputweight, options.f, options.df)     
    else
        error("Model " * options.modeltype * " is not supported")
    end        
    
    options.featurenorm = initweight

    return model, options     
end
    
function LossFunction(model, data, options::ModelOptions)
    loss = 0.0
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, options.batchsize) 
        lossj = LossFunction(model, options, X, Y) 
        loss = loss + lossj     
    end
    return loss 
end

function LossFunctionEnzyme(model, data, options::ModelOptions)
    loss = 0.0
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, options.batchsize) 
        lossj = LossFunctionEnzyme(model, options, X, Y) 
        loss = loss + lossj     
    end
    return loss 
end

function GetDescriptors(data, options::ModelOptions)    
    D = Array{Any}(nothing, length(data.fileindex))
    for i = 1:length(data.fileindex)
        X, Y = LoadData(data, i, options.batchsize)         
        D[i] = X[1]
    end    
    return D
end

function UpdateModel(model, data, options::ModelOptions, ipoch)

    epsilon = 0e-2 
    alpha = options.η*(1 + epsilon*(ipoch-1))     
    for j = 1:length(data.fileindex)
        X, Y = LoadData(data, j, options.batchsize) 
        gradj = LossGradient(model, options, X, Y)  
        # model.convweight = model.convweight - alpha[1]*gradj[1]
        # model.convbias = model.convbias - alpha[2]*gradj[2]
        # model.outputweight = model.outputweight - alpha[3]*gradj[3]
        if sum(isnan.(gradj[1])) > 0             
            @warn "Nan detected"      
        else
            model.convweight = model.convweight - alpha[1]*gradj[1]
        end
        if sum(isnan.(gradj[2])) > 0 
            @warn "Nan detected"      
        else
            model.convbias = model.convbias - alpha[2]*gradj[2]
        end
        if sum(isnan.(gradj[3])) > 0 
            @warn "Nan detected"      
        else
            model.outputweight = model.outputweight - alpha[3]*gradj[3]
        end
        if options.modeltype == "diffnn"     
            if sum(isnan.(gradj[4])) > 0 
                @warn "Nan detected"      
            else
                model.activationweight = model.activationweight - alpha[4]*gradj[4]
            end 
            if sum(isnan.(gradj[5])) > 0 
                @warn "Nan detected"      
            else
                model.activationscale = model.activationscale - alpha[5]*gradj[5]
            end 
        end        
        # model.convweight = model.convweight - alpha[1]*gradj[1]
        # model.convbias = model.convbias - alpha[2]*gradj[2]
        # model.outputweight = model.outputweight - alpha[3]*gradj[3]
        # if options.modeltype == "diffnn"      
        #     model.activationweight = model.activationweight - alpha[4]*gradj[4]
        #     model.activationscale = model.activationscale - alpha[5]*gradj[5]
        # end
    end                     
    return model
end

function TrainModel(model, data, options::ModelOptions)
    loss = LossFunction(model, data, options)
    display("Loss value: $loss")
    for i = 1:options.epochs
        model = UpdateModel(model, data, options, i)
        loss = LossFunction(model, data, options)
        maeenergy, maeforce = MAE(model, data, options)
        display("Loss value: $loss")
    end
    return model 
end

using Potential 

function TrainModel(traindata, descriptors, Doptions, modelOpts::ModelOptions, potential=nothing) 

    filename = modelOpts.filename;
    modelOpts.normalizeenergy = Doptions.normalizeenergy
    
    if (modelOpts.modeltype == "dnn") || (modelOpts.modeltype == "diffnn")
        # Generate training data for DNN model                 
        data = Potential.DNNdata(traindata, descriptors, Doptions, filename, modelOpts.batchsize, potential) 
    elseif modelOpts.modeltype == "gnn"
        data = Potential.GNNdata(traindata, descriptors, Doptions, filename, modelOpts.batchsize, potential) 
    else
        error("Model " * options.modeltype * " is not supported")
    end    

    # Build a DML model 
    model, modelOpts = BuildModel(data, modelOpts)

    # display(model.convweight)
    # display(model.convbias)
    # display(model.outputweight)
    # display(model.activationweight)
    # display(model.activationscale)

    # calculate errors 
    maeenergy, maeforce = MAE(model, data, modelOpts)

    # Train the DML model
    model = TrainModel(model, data, modelOpts)

    return model, data, modelOpts  
end

function MAE(model, data, options::ModelOptions)
    maeenergy = 0.0
    maeforce = 0.0
    k = 0;
    n = 0;
    for i = 1:length(data.fileindex)
        x, y = LoadData(data, i, options.batchsize) 
        if (options.modeltype=="dnn") || (options.modeltype=="diffnn")
            for j=1:size(x[1],2)
                D = x[1][:,j]
                dDdR = x[2][:,:,j]
                energy, forces, dedD = model(D, dDdR);                            
                enorm = 3.0/length(forces)     
                if options.normalizeenergy == true
                    enorm = 1.0
                end
                maeenergy = maeenergy + abs(energy-y[1][j])*enorm                      
                tm = y[2][:,:,j] 
                maeforce = maeforce + sum(abs.(forces[:] - tm[:]))
                k = k + 1;
                n = n + length(forces)
            end    
        elseif options.modeltype=="gnn"
            for j=1:size(x[1],3)
                D = x[1][:,:,j]
                dDdR = x[2][:,:,:,j]
                A = x[3][:,:,j]
                P = x[4][:,:,j]    
                energy, forces, dedD = model(D, A, P, dDdR);         
                enorm = 3.0/length(forces)     
                if options.normalizeenergy == true
                    enorm = 1.0
                end       
                maeenergy = maeenergy + abs(energy-y[1][j])*enorm
                tm = y[2][:,:,j] 
                maeforce = maeforce + sum(abs.(forces[:] - tm[:]))
                k = k + 1;
                n = n + length(forces)
            end    
        end        
    end
    maeenergy = maeenergy/k
    maeforce = maeforce/n 
    print("Energy MAE Error: $maeenergy \n")
    print("Force MAE Error: $maeforce \n")
    return maeenergy, maeforce, k 
end


# function energy(model, D)
#     return  model(D)                
# end
# function energyforces(model, D, dDdR)
#     return  model(D, dDdR)                
# end
