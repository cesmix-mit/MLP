function lossenergy(Epredict, Edata)
    loss = (Epredict - Edata)^2
    return loss     
end

function lossenergy(Epredict, Edata, forcelength)
    loss = (Epredict*(3.0/forcelength) - Edata*(3.0/forcelength))^2
    return loss     
end

function lossse(Epredict, Fpredict, Edata, Fdata, β, forcelength)
    loss = β[1]*(Epredict - Edata)^2
    for i = 1:forcelength
        loss = loss + β[2]*(Fpredict[i] - Fdata[i])^2 
    end    
    return loss     
end

function lossmse(Epredict, Fpredict, Edata, Fdata, β, forcelength)
    α = β[2]/forcelength
    loss = β[1]*(Epredict*(3.0/forcelength) - Edata*(3.0/forcelength))^2
    for i = 1:forcelength
        loss = loss + α*(Fpredict[i] - Fdata[i])^2 
    end    
    return loss     
end

function lossmsenorm(Epredict, Fpredict, Edata, Fdata, β, forcelength)
    α = β[2]/forcelength
    loss = β[1]*(Epredict - Edata)^2
    for i = 1:forcelength
        loss = loss + α*(Fpredict[i] - Fdata[i])^2 
    end    
    return loss     
end

function lossae(Epredict, Fpredict, Edata, Fdata, β, forcelength)
    loss = β[1]*abs(Epredict*(3.0/forcelength) - Edata*(3.0/forcelength))    
    for i = 1:forcelength
        loss = loss + β[2]*abs(Fpredict[i] - Fdata[i])
    end
    return loss 
end

function lossmae(Epredict, Fpredict, Edata, Fdata, β, forcelength)
    α = β[2]/(forcelength)
    loss = β[1]*abs(Epredict*(3.0/forcelength) - Edata*(3.0/forcelength))    
    for i = 1:forcelength
        loss = loss + α*abs(Fpredict[i] - Fdata[i])
    end
    return loss 
end

# function lossmse(Yhat, Y, β)
#     Ehat = Yhat[1] # Energy prediction from a DML model 
#     Fhat = Yhat[2] # Force prediction from a DML model 
#     E = Y[1]       # Energy prediction from DFT calculation  
#     F = Y[2]       # Force prediction from DFT calculation  
#     loss = (1.0-β)*sum((Ehat[:] - E[:]).^2) + β*sum((Fhat[:] - F[:]).^2) 
#     return loss 
# end

