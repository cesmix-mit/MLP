using Preprocessing, Potential


function linearsolve(Ae, Af, As, be, bf, bs)

    cefs = zeros(0)
    if (length(Ae)>0) & (length(Af)>0) & (length(As)>0) & (length(be)>0) & (length(bf)>0) & (length(bs)>0) 
        cefs = cat(cat(Ae, Af, dims=1), As, dims=1)\[be[:]; bf[:]; bs[:]]
    end

    cef = zeros(0)
    if (length(Ae)>0) & (length(Af)>0) & (length(be)>0) & (length(bf)>0) 
        cef = cat(Ae, Af, dims=1)\[be[:]; bf[:]]                
    end

    ce = zeros(0)
    if (length(Ae)>0) & (length(be)>0)
        ce = Ae\be
    end

    cf = zeros(0)
    if (length(Af)>0) & (length(bf)>0)
        cf = Af\bf                
    end

    return ce, cf, cef, cefs
end

function lsqsolve(Ae, Af, As, be, bf, bs)

    cefs = zeros(0)
    cef = zeros(0)
    ce = zeros(0)
    cf = zeros(0)

    sumAe = sum(abs.(Ae[:])) > 1e-15
    sumAf = sum(abs.(Af[:])) > 1e-15
    sumAs = sum(abs.(As[:])) > 1e-15
    sumbe = sum(abs.(be[:])) > 1e-15
    sumbf = sum(abs.(bf[:])) > 1e-15
    sumbs = sum(abs.(bs[:])) > 1e-15

    if (sumAe) & (sumAf) & (sumAs) & (sumbe) & (sumbf) & (sumbs) 
        A = Ae+Af+As
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end
        cefs = (A)\(be+bf+bs)
    end
    if (sumAe) & (sumAf) & (sumbe) & (sumbf) 
        A = Ae+Af
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end
        cef = A\(be+bf)
    end
    if (sumAe) & (sumbe) 
        A = 1.0*Ae
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end
        ce = A\be
    end
    if (sumAf) & (sumbf) 
        A = 1.0*Af
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end        
        if abs(A[1,1])<=1e-10
            A[1,1] = 1e-10
        end
        cf = A\bf
    end

    return ce, cf, cef, cefs
end

function insertarray(ce, cei, i)
    if (length(cei)>0)
        ce[:,i] = cei
    end
    return ce
end

function catarray(Ae, be, Aei, bei)
    if (length(Aei) > 0) & (length(bei) > 0)
        Ae = cat(Ae, Aei, dims=1)
        be = [be; bei]
    end
    return Ae, be
end

function calmae(Aei, bei, cei)
    mae = -1.0
    if (length(bei) > 0) & (length(cei)>0) & (length(Aei)>0)       
        mae = sum(abs.(bei - Aei*cei))/length(bei)       
    end
    return mae
end

function calerrors(Aei, bei, cei)
    mae = -1.0
    rmse = - 1.0
    rsq = -1.0
    ssr = -1.0
    avg = -1.0
    if (length(bei) > 0) & (length(cei)>0) & (length(Aei)>0)       
        res = bei - Aei*cei
        mae = sum(abs.(res))/length(bei)       
        ssr = sum(res.^2)
        mse = ssr /length(bei);
        rmse = sqrt(mse) 
        avg = sum(bei)/length(bei);
        rsq = 1.0 - ssr/sum((bei .- avg).^2);
    end
    return mae, rmse, rsq, ssr, avg
end

function applyweight(Aei, bei, we)
    if (length(bei) > 0) & (length(Aei) > 0) & (length(we) > 0) 
        Aei = we[1]*Aei
        bei = we[1]*bei
    end
    return Aei, bei 
end

function applytranspose(Aei, bei, m, n)    
    Cei = zeros(m,n)
    if (length(bei) > 0) & (length(Aei) > 0) 
        Cei = Aei'*bei        
    end
    return Cei
end

function linearfit(data, descriptors, potential, Doptions, optim::OptimStruct)
        
    data = Preprocessing.deletenothing(data)
    descriptors = Preprocessing.deletenothing(descriptors)
    potential = Preprocessing.deletenothing(potential)
    
    lossfunc, method, eta, kappa, weightouter, etaspace, kappalist = getoptim(optim)

    pbc = Doptions.pbc
    normalizeenergy = Doptions.normalizeenergy
    normalizestress = Doptions.normalizestress
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c
    
    n = length(data)
    emae = -ones(Float64,n,4)
    fmae = -ones(Float64,n,4)
    smae = -ones(Float64,n,4)
    ce = 1.0
    cf = 1.0
    cef = 1.0
    cefs = 1.0
    Ae = 1.0
    Af = 1.0
    As = 1.0
    be = 1.0
    bf = 1.0
    bs = 1.0
    De = 1.0
    Df = 1.0
    Ds = 1.0
    de = 1.0
    df = 1.0
    ds = 1.0
    Le = 1.0
    Lf = 1.0
    Ls = 1.0
    he = 1.0
    hf = 1.0
    hs = 1.0
    for i = 1:n
        # training data set
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        training = Array(1:config.nconfigs)      

        # display(config.we)
        # display(config.wf)
        # display(config.ws)
        # error("here")
    
        Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.globaldescriptors(config, training, eta, kappa, 
                normalizeenergy, normalizestress, descriptors, potential)

        m = size(Aei,2) # total number of descriptors
        if i == 1            
            ce = ce*zeros(Float64,m,n)
            cf = cf*zeros(Float64,m,n)
            cef = cef*zeros(Float64,m,n)
            cefs = cefs*zeros(Float64,m,n)
            De = De*zeros(m,m,n)
            Df = Df*zeros(m,m,n)
            Ds = Ds*zeros(m,m,n)
            de = de*zeros(m,n)
            df = df*zeros(m,n)
            ds = ds*zeros(m,n)
            Le = Le*zeros(m,m)
            Lf = Lf*zeros(m,m)
            Ls = Ls*zeros(m,m)
            he = he*zeros(m)
            hf = hf*zeros(m)
            hs = hs*zeros(m)
        end
        
        De[:,:,i] = applytranspose(Aei, Aei, m, m)    
        Df[:,:,i] = applytranspose(Afi, Afi, m, m)    
        Ds[:,:,i] = applytranspose(Asi, Asi, m, m)    
        de[:,i] = applytranspose(Aei, bei, m, 1)    
        df[:,i] = applytranspose(Afi, bfi, m, 1)    
        ds[:,i] = applytranspose(Asi, bsi, m, 1)    

        Aei, bei = applyweight(Aei, bei, config.we)
        Afi, bfi = applyweight(Afi, bfi, config.wf)
        Asi, bsi = applyweight(Asi, bsi, config.ws)
                
        Le += applytranspose(Aei, Aei, m, m)    
        Lf += applytranspose(Afi, Afi, m, m)    
        Ls += applytranspose(Asi, Asi, m, m)    
        he += applytranspose(Aei, bei, m, 1)    
        hf += applytranspose(Afi, bfi, m, 1)    
        hs += applytranspose(Asi, bsi, m, 1)            
        
        if length(data)==1
            if optim.method == "lsq"
                cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
            else
                cei, cfi, cefi, cefsi = linearsolve(Aei, Afi, Asi, bei, bfi, bsi)                        
            end            
            ce = insertarray(ce, cei, i)
            cf = insertarray(cf, cfi, i)
            cef = insertarray(cef, cefi, i)
            cefs = insertarray(cefs, cefsi, i)

            emae[i,1] = calmae(Aei, bei, cei)
            emae[i,2] = calmae(Aei, bei, cfi)
            emae[i,3] = calmae(Aei, bei, cefi)
            emae[i,4] = calmae(Aei, bei, cefsi)
            fmae[i,1] = calmae(Afi, bfi, cei)
            fmae[i,2] = calmae(Afi, bfi, cfi)
            fmae[i,3] = calmae(Afi, bfi, cefi)
            fmae[i,4] = calmae(Afi, bfi, cefsi)
            smae[i,1] = calmae(Asi, bsi, cei)
            smae[i,2] = calmae(Asi, bsi, cfi)
            smae[i,3] = calmae(Asi, bsi, cefi)
            smae[i,4] = calmae(Asi, bsi, cefsi)

            if lossfunc == "energy" 
                coeff = ce
            elseif lossfunc == "force" 
                coeff = cf
            elseif lossfunc == "energyforce" 
                coeff = cef
            elseif lossfunc == "energyforcestress" 
                coeff = cefs
            end    
                    
            return coeff, ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
        end

        if (i == 1) | (optim.method == "lsq")
            Ae = 1.0*Aei
            Af = 1.0*Afi
            As = 1.0*Asi
            be = 1.0*bei
            bf = 1.0*bfi
            bs = 1.0*bsi
        else         
            Ae, be = catarray(Ae, be, Aei, bei)
            Af, bf = catarray(Af, bf, Afi, bfi)
            As, bs = catarray(As, bs, Asi, bsi)                        
        end                        

        if (optim.method == "lsq")
            #display(Le)
            cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
        else
            cei, cfi, cefi, cefsi = linearsolve(Ae, Af, As, be, bf, bs)
        end

        ce = insertarray(ce, cei, i)
        cf = insertarray(cf, cfi, i)
        cef = insertarray(cef, cefi, i)
        cefs = insertarray(cefs, cefsi, i)

        emae[i,1] = calmae(Ae, be, cei)
        emae[i,2] = calmae(Ae, be, cfi)
        emae[i,3] = calmae(Ae, be, cefi)
        emae[i,4] = calmae(Ae, be, cefsi)
        fmae[i,1] = calmae(Af, bf, cei)
        fmae[i,2] = calmae(Af, bf, cfi)
        fmae[i,3] = calmae(Af, bf, cefi)
        fmae[i,4] = calmae(Af, bf, cefsi)
        smae[i,1] = calmae(As, bs, cei)
        smae[i,2] = calmae(As, bs, cfi)
        smae[i,3] = calmae(As, bs, cefi)
        smae[i,4] = calmae(As, bs, cefsi)
    end

    if lossfunc == "energy" 
        coeff = ce[:,end]
    elseif lossfunc == "force" 
        coeff = cf[:,end]
    elseif lossfunc == "energyforce" 
        coeff = cef[:,end]
    elseif lossfunc == "energyforcestress" 
        coeff = cefs[:,end]
    end    
    
    return coeff, ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
end

function validate(data, descriptors, potential, Doptions, optim::OptimStruct, coeff)    
    
    data = Preprocessing.deletenothing(data)
    descriptors = Preprocessing.deletenothing(descriptors)
    potential = Preprocessing.deletenothing(potential)

    lossfunc, method, eta, kappa, weightouter, etaspace, kappalist = getoptim(optim)
    
    pbc = Doptions.pbc
    normalizeenergy = Doptions.normalizeenergy
    normalizestress = Doptions.normalizestress
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c
    
    n = length(data)
    emae = -ones(Float64,n)
    fmae = -ones(Float64,n)
    smae = -ones(Float64,n)
    ermse = -ones(Float64,n)
    frmse = -ones(Float64,n)
    srmse = -ones(Float64,n)
    ersq = -ones(Float64,n)
    frsq = -ones(Float64,n)
    srsq = -ones(Float64,n)
    essr = -ones(Float64,n)
    fssr = -ones(Float64,n)
    sssr = -ones(Float64,n)
    eavg = -ones(Float64,n)
    favg = -ones(Float64,n)
    savg = -ones(Float64,n)
    szbe = zeros(Int64, n)
    szbf = zeros(Int64, n)
    szbs = zeros(Int64, n)
    be = []
    bf = []
    bs = []
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        valid = Array(1:config.nconfigs)      

        Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.globaldescriptors(config, valid, eta, kappa, 
                        normalizeenergy, normalizestress, descriptors, potential)

        szbe[i] = length(bei)        
        szbf[i] = length(bfi)        
        szbs[i] = length(bsi)        

        emae[i], ermse[i], ersq[i], essr[i], eavg[i] = calerrors(Aei, bei, coeff)
        fmae[i], frmse[i], frsq[i], fssr[i], favg[i] = calerrors(Afi, bfi, coeff)
        smae[i], srmse[i], srsq[i], sssr[i], savg[i] = calerrors(Asi, bsi, coeff)

        # display(sum(abs.(bei))/length(bei))
        # display(emae[i])
        # display(sum(abs.(bfi))/length(bfi))
        # display(fmae[i])

        if szbe[i] > 0
            be = [be; bei]
        end 
        if szbf[i] > 0
            bf = [bf; bfi]
        end 
        if szbs[i] > 0
            bs = [bs; bsi]
        end 
    end

    energyerrors = -ones((n+1),3)    
    forceerrors = -ones((n+1),3)    
    stresserrors = -ones((n+1),3)    

    if sum(szbe) > 0 
        ssr = sum(essr)
        avg = sum(eavg.*szbe)/sum(szbe)          
        rsq = 1.0 - ssr/sum((be .- avg).^2);
        energyerrors[1,1] = sum(emae.*szbe)/sum(szbe)  # MAE
        energyerrors[1,2] = sqrt(sum((ermse.^2).*szbe)/sum(szbe)) # RMSE 
        energyerrors[1,3] = rsq
        energyerrors[2:end,:] = [emae ermse ersq]        
    end

    if sum(szbf) > 0 
        ssr = sum(fssr)
        avg = sum(favg.*szbf)/sum(szbf)          
        rsq = 1.0 - ssr/sum((bf .- avg).^2);
        forceerrors[1,1] = sum(fmae.*szbf)/sum(szbf)
        forceerrors[1,2] =  sqrt(sum((frmse.^2).*szbf)/sum(szbf))    
        forceerrors[1,3] = rsq 
        forceerrors[2:end,:] = [fmae frmse frsq]   
    end    
    
    if sum(szbs) > 0 
        ssr = sum(sssr)
        avg = sum(savg.*szbs)/sum(szbs)          
        rsq = 1.0 - ssr/sum((bs .- avg).^2);
        stresserrors[1,1] = sum(smae.*szbs)/sum(szbs)
        stresserrors[1,2] = sqrt(sum((srmse.^2).*szbs)/sum(szbs))        
        stresserrors[1,3] = rsq 
        stresserrors[2:end,:] = [smae srmse srsq]   
    end

    return energyerrors, forceerrors, stresserrors
end

function setcutoff(descriptors, eta)
    if length(eta)>0
        for i = 1:length(descriptors)    
            if typeof(descriptors[i]) == Potential.PotentialStruct
                descriptors[i].rcut = eta[1] 
            elseif typeof(descriptors[i]) == Potential.SnaStruct
                descriptors[i].rcutfac = eta[1]
            elseif typeof(descriptors[i]) == Potential.SNAPparams
                descriptors[i].rcutfac = eta[1]                
                #descriptors[i] = initsnap(descriptors[i])    
            elseif descriptors[i].name == "ACE" 
                descriptors[i].rcut = eta[1]
                if length(eta)>1
                    descriptors[i].rin = eta[2]          
                end       
            elseif descriptors[i].name == "POD" 
                descriptors[i].rcut = eta[1]          
                if length(eta)>1
                    descriptors[i].rin = eta[2]          
                end       
                descriptors[i] = initPOD(descriptors[i])
            end            
        end
    end
    return descriptors
end

function maespace(traindata, validdata, descriptors, potential, Doptions, optim::OptimStruct)

    # generate the eta space grid
    etaspace = optim.etaspace
    etamin = etaspace[:,1]
    etamax = etaspace[:,2]
    N = Int64.(etaspace[:,3])
    eta = tensornodes(etamin, etamax, N)
        
    pbc = Doptions.pbc
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c

    # read data 
    trainconfig = Array{Any}(undef, length(traindata))
    for i = 1:length(traindata)
        trainconfig[i], indices = Preprocessing.readconfigdata(traindata[i], pbc, a, b, c)                  
    end
    validconfig = Array{Any}(undef, length(validdata))
    for i = 1:length(validdata)
        validconfig[i], indices = Preprocessing.readconfigdata(validdata[i], pbc, a, b, c)                  
    end

    neta = size(eta,1)
    meta = size(eta,2)
    emae = -ones(neta)
    fmae = -ones(neta)
    smae = -ones(neta)
    for n = 1:neta
        optim.eta[1:meta] = eta[n,:]
        descriptors = setcutoff(descriptors, optim.eta)      
        @time begin
        coeff,~ = linearfit(trainconfig, descriptors, potential, Doptions, optim)
        if n==1
            display(size(coeff))
        end
        eerr, ferr, serr = validate(validconfig, descriptors, potential, Doptions, optim, coeff)    
        end
        emae[n] = eerr[1]
        fmae[n] = ferr[1]
        smae[n] = serr[1]                
        a1 = num2string(emae[n], 16); 
        a2 = num2string(fmae[n], 16); 
        a3 = num2string(smae[n], 16);
        mystr = "$n  "
        for j = 1:length(eta[n,:])
            a = num2string(eta[n,j], 16)
            mystr = mystr * a * "  "
        end
        print(mystr * a1 * "  " * a2 * "  " * a3 * " \n")
    end

    return eta, emae, fmae, smae
end

function gradientdescent(c, yspace, y0=nothing)

    dim = size(yspace,1)
    ymin = yspace[:,1]
    ymax = yspace[:,2]
    N = Int64.(yspace[:,3])
    
    if y0 === nothing
        M = 100;
        y0 = rand(M, dim)
        for i = 1:dim
            y0[:,i] = ymin[i] .+ y0[:,i]*(ymax[i] - ymin[i])
        end
    end
    
    # estimate maximum stepsize
    f, df = mpolyder(y0, c, ymin, ymax, N)
    y = sum(y0.*y0, dims=2) 
    gamma = sum(df.*df, dims=2) 
    alphamax = minimum(sqrt.(y./gamma))
    alphamax = max(alphamax, 10.0)

    beta = 0.8;    
    y = 1.0*y0
    M = size(y,1) 
    y1 = 1.0*y
    r1 = zeros(M)
    maxstep = zeros(M)    
    alpha = alphamax*ones(M)
    a = zeros(dim)
    iter = 0
    while (true)
        iter = iter + 1
        f, df = mpolyder(y, c, ymin, ymax, N)        
        gamma = 0.5*sum(df.*df, dims=2) 
        
        if (minimum(gamma[:]) <= 1e-12) || (iter >= 200)
            fmin, imin = findmin(f)
            sol = y[imin,:]
            return sol, fmin, iter, y, f, df, y0  
        end

        for m = 1:M
            for i = 1:dim
                a1 = (y[m,i]-ymax[i])/df[m,i]
                a2 = (y[m,i]-ymin[i])/df[m,i]    
                a[i] = (abs(df[m,i])<=1e-12) ? alphamax : max(a1, a2)    
            end
            maxstep[m] = min(minimum(a), alphamax)
            alpha[m] = min(alpha[m], maxstep[m])
            y1[m,:] = y[m,:] - alpha[m]*df[m,:]    
            r1[m] = f[m] - alpha[m]*gamma[m]
        end
                        
        f1 = mpolyeval(y1, c, ymin, ymax, N)
        for m = 1:M
            fm = f1[m]
            rm = r1[m]              
            incm = 0              
            while (fm > rm)
                alpha[m] = beta*alpha[m]
                rm = f[m] - alpha[m]*gamma[m]
                y1[m,:] = y[m,:] - alpha[m]*df[m,:]                    
                tm = mpolyeval(reshape(y1[m,:],(1,dim)), c, ymin, ymax, N)                
                fm = tm[1]
                incm += 1 
                if incm > 10 
                   break; 
                end
            end            
            y[m,:] = y1[m,:]                         
            alpha[m] = (1.0/beta)*alpha[m]
        end        
    end        
end

function optimize(traindata, testdata, descriptors, potential, Doptions, optim::OptimStruct)

    traindata = Preprocessing.deletenothing(traindata)
    testdata = Preprocessing.deletenothing(testdata)
    descriptors = Preprocessing.deletenothing(descriptors)
    potential = Preprocessing.deletenothing(potential)

    # compute errors on the eta space grid
    etapts, emae, fmae, smae = maespace(traindata, testdata, descriptors, potential, Doptions, optim)    

    lossfunc = optim.lossfunc
    weightouter = optim.weightouter
    etaspace = optim.etaspace

    # compute the loss function on the eta space grid
    if lossfunc == "energy" 
        f = weightouter[1]*emae
    elseif lossfunc == "force" 
        f = weightouter[2]*emae
    elseif lossfunc == "energyforce" 
        f = weightouter[1]*emae + weightouter[2]*fmae
    elseif lossfunc == "energyforcestress" 
        f = weightouter[1]*emae + weightouter[2]*fmae + weightouter[3]*smae
    else
        error("lossfunc is not valid. It must be energy, force, energyforce, or energyforcestress")
    end    

    # interpolate the loss function on the eta space grid
    pc = mpolyinterp(etapts, f, etaspace)

    # optimize the loss function on the eta space
    eta, fmin, iter  = gradientdescent(pc, etaspace)
    meta = length(eta)

    if fmin > minimum(f[:])
        fmn,ind = findmin(f[:])
        eta = etapts[ind,:]
    end

    optim.eta[1:meta] = eta

    # set cut-off radius
    descriptors = setcutoff(descriptors, optim.eta)        
    
    # compute coefficient 
    coeff,~ = linearfit(traindata, descriptors, potential, Doptions, optim)

    return eta, coeff, fmin, iter, pc, f, etapts, emae, fmae, smae
end


# function calmean(emae, szbe)
#     n,p,m = size(emae)
#     em = -ones(p,m)       
#     if (sum(szbe) > 0) & (sum(emae[:]) > 0)
#         for k = 1:p
#             for j = 1:m
#                 em[k,j] = 0.0
#                 for i = 1:n            
#                     em[k,j] = em[k,j] + emae[i,k,j]*szbe[i]
#                 end
#                 em[k,j] =  em[k,j]/sum(szbe)
#             end
#         end
#     end
#     return em
# end

# function linearfit(data, emdescriptors, mldescriptors, eta, kappa, pbc=nothing, a=nothing, b=nothing, c=nothing, method=nothing, normalizeenergy=nothing)
    
#     if normalizeenergy === nothing
#         normalizeenergy = 0 # divide energy by number of atoms
#     end
#     if method === nothing
#         method = "lsq"
#     end

#     n = length(data)
#     emae = -ones(Float64,n,4)
#     fmae = -ones(Float64,n,4)
#     smae = -ones(Float64,n,4)
#     ce = 1.0
#     cf = 1.0
#     cef = 1.0
#     cefs = 1.0
#     Ae = 1.0
#     Af = 1.0
#     As = 1.0
#     be = 1.0
#     bf = 1.0
#     bs = 1.0
#     De = 1.0
#     Df = 1.0
#     Ds = 1.0
#     de = 1.0
#     df = 1.0
#     ds = 1.0
#     Le = 1.0
#     Lf = 1.0
#     Ls = 1.0
#     he = 1.0
#     hf = 1.0
#     hs = 1.0
#     for i = 1:n
#         # training data set
#         if typeof(data[i]) == Preprocessing.DataStruct
#             config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
#         else
#             config = data[i]
#         end
#         training = Array(1:config.nconfigs)      

#         Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.emldescriptors(config, training, eta, kappa, normalizeenergy, emdescriptors, mldescriptors)
        
#         m = size(Aei,2) # total number of descriptors
#         if i == 1            
#             ce = ce*zeros(Float64,m,n)
#             cf = cf*zeros(Float64,m,n)
#             cef = cef*zeros(Float64,m,n)
#             cefs = cefs*zeros(Float64,m,n)
#             De = De*zeros(m,m,n)
#             Df = Df*zeros(m,m,n)
#             Ds = Ds*zeros(m,m,n)
#             de = de*zeros(m,n)
#             df = df*zeros(m,n)
#             ds = ds*zeros(m,n)
#             Le = Le*zeros(m,m)
#             Lf = Lf*zeros(m,m)
#             Ls = Ls*zeros(m,m)
#             he = he*zeros(m)
#             hf = hf*zeros(m)
#             hs = hs*zeros(m)
#         end
        
#         De[:,:,i] = applytranspose(Aei, Aei, m, m)    
#         Df[:,:,i] = applytranspose(Afi, Afi, m, m)    
#         Ds[:,:,i] = applytranspose(Asi, Asi, m, m)    
#         de[:,i] = applytranspose(Aei, bei, m, 1)    
#         df[:,i] = applytranspose(Afi, bfi, m, 1)    
#         ds[:,i] = applytranspose(Asi, bsi, m, 1)    

#         Aei, bei = applyweight(Aei, bei, config.we)
#         Afi, bfi = applyweight(Afi, bfi, config.wf)
#         Asi, bsi = applyweight(Asi, bsi, config.ws)
        
#         Le += applytranspose(Aei, Aei, m, m)    
#         Lf += applytranspose(Afi, Afi, m, m)    
#         Ls += applytranspose(Asi, Asi, m, m)    
#         he += applytranspose(Aei, bei, m, 1)    
#         hf += applytranspose(Afi, bfi, m, 1)    
#         hs += applytranspose(Asi, bsi, m, 1)    

#         if length(data)==1
#             if method == "lsq"
#                 cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
#             else
#                 cei, cfi, cefi, cefsi = linearsolve(Aei, Afi, Asi, bei, bfi, bsi)                        
#             end            
#             ce = insertarray(ce, cei, i)
#             cf = insertarray(cf, cfi, i)
#             cef = insertarray(cef, cefi, i)
#             cefs = insertarray(cefs, cefsi, i)

#             emae[i,1] = calmae(Aei, bei, cei)
#             emae[i,2] = calmae(Aei, bei, cfi)
#             emae[i,3] = calmae(Aei, bei, cefi)
#             emae[i,4] = calmae(Aei, bei, cefsi)
#             fmae[i,1] = calmae(Afi, bfi, cei)
#             fmae[i,2] = calmae(Afi, bfi, cfi)
#             fmae[i,3] = calmae(Afi, bfi, cefi)
#             fmae[i,4] = calmae(Afi, bfi, cefsi)
#             smae[i,1] = calmae(Asi, bsi, cei)
#             smae[i,2] = calmae(Asi, bsi, cfi)
#             smae[i,3] = calmae(Asi, bsi, cefi)
#             smae[i,4] = calmae(Asi, bsi, cefsi)
            
#             return ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
#         end

#         if (i == 1) | (method == "lsq")
#             Ae = Aei
#             Af = Afi
#             As = Asi
#             be = bei
#             bf = bfi
#             bs = bsi
#         else         
#             Ae, be = catarray(Ae, be, Aei, bei)
#             Af, bf = catarray(Af, bf, Afi, bfi)
#             As, bs = catarray(As, bs, Asi, bsi)
#         end
        
#         if (method == "lsq")
#             cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
#         else
#             cei, cfi, cefi, cefsi = linearsolve(Ae, Af, As, be, bf, bs)
#         end

#         ce = insertarray(ce, cei, i)
#         cf = insertarray(cf, cfi, i)
#         cef = insertarray(cef, cefi, i)
#         cefs = insertarray(cefs, cefsi, i)

#         emae[i,1] = calmae(Ae, be, cei)
#         emae[i,2] = calmae(Ae, be, cfi)
#         emae[i,3] = calmae(Ae, be, cefi)
#         emae[i,4] = calmae(Ae, be, cefsi)
#         fmae[i,1] = calmae(Af, bf, cei)
#         fmae[i,2] = calmae(Af, bf, cfi)
#         fmae[i,3] = calmae(Af, bf, cefi)
#         fmae[i,4] = calmae(Af, bf, cefsi)
#         smae[i,1] = calmae(As, bs, cei)
#         smae[i,2] = calmae(As, bs, cfi)
#         smae[i,3] = calmae(As, bs, cefi)
#         smae[i,4] = calmae(As, bs, cefsi)
#     end

#     return ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
# end

# function validation(data, emdescriptors, mldescriptors, ce, cf, cef, cefs, eta, kappa, pbc=nothing, a=nothing, b=nothing, c=nothing)    

#     normalizeenergy = 0
#     m = size(ce,2)
#     n = length(data)
#     emae = -ones(Float64,n,4,m)
#     fmae = -ones(Float64,n,4,m)
#     smae = -ones(Float64,n,4,m)
#     szbe = zeros(Int64, n)
#     szbf = zeros(Int64, n)
#     szbs = zeros(Int64, n)
#     for i = 1:n
#         if typeof(data[i]) == Preprocessing.DataStruct
#             config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
#         else
#             config = data[i]
#         end
#         valid = Array(1:config.nconfigs)      

#         Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.emldescriptors(config, valid, eta, kappa, normalizeenergy, emdescriptors, mldescriptors)
#         szbe[i] = length(bei)        
#         szbf[i] = length(bfi)        
#         szbs[i] = length(bsi)        

#         for j = 1:m
#             cei = ce[:,j]
#             cfi = cf[:,j]
#             cefi = cef[:,j]
#             cefsi = cefs[:,j] 
#             emae[i,1,j] = calmae(Aei, bei, cei)
#             emae[i,2,j] = calmae(Aei, bei, cfi)
#             emae[i,3,j] = calmae(Aei, bei, cefi)
#             emae[i,4,j] = calmae(Aei, bei, cefsi)
#             fmae[i,1,j] = calmae(Afi, bfi, cei)
#             fmae[i,2,j] = calmae(Afi, bfi, cfi)
#             fmae[i,3,j] = calmae(Afi, bfi, cefi)
#             fmae[i,4,j] = calmae(Afi, bfi, cefsi)
#             smae[i,1,j] = calmae(Asi, bsi, cei)
#             smae[i,2,j] = calmae(Asi, bsi, cfi)
#             smae[i,3,j] = calmae(Asi, bsi, cefi)
#             smae[i,4,j] = calmae(Asi, bsi, cefsi)
#         end
#     end

#     eme = calmean(emae, szbe)
#     fme = calmean(fmae, szbf)
#     sme = calmean(smae, szbs)

#     return eme, fme, sme, emae, fmae, smae, szbe, szbf, szbs
# end

# function linearfit2(data, descriptors, potential, optim, pbc=nothing, a=nothing, b=nothing, c=nothing)
        
#     data = Preprocessing.deletenothing(data)
#     descriptors = Preprocessing.deletenothing(descriptors)
#     potential = Preprocessing.deletenothing(potential)
    
#     if pbc === nothing        
#         pbc = [1; 1; 1] # periodic boundary conditions
#     end

#     lossfunc, method, normalizeenergy, normalizestress, eta, kappa, 
#         weightouter, etaspace, kappalist = getoptim(optim)

#     n = length(data)
#     emae = -ones(Float64,n,4)
#     fmae = -ones(Float64,n,4)
#     smae = -ones(Float64,n,4)
#     ce = 1.0
#     cf = 1.0
#     cef = 1.0
#     cefs = 1.0
#     Ae = 1.0
#     Af = 1.0
#     As = 1.0
#     be = 1.0
#     bf = 1.0
#     bs = 1.0
#     De = 1.0
#     Df = 1.0
#     Ds = 1.0
#     de = 1.0
#     df = 1.0
#     ds = 1.0
#     Le = 1.0
#     Lf = 1.0
#     Ls = 1.0
#     he = 1.0
#     hf = 1.0
#     hs = 1.0
#     for i = 1:n
#         # training data set
#         if typeof(data[i]) == Preprocessing.DataStruct
#             config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
#         else
#             config = data[i]
#         end
#         training = Array(1:config.nconfigs)      

#         Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.globaldescriptors(config, training, eta, kappa, 
#                 normalizeenergy, normalizestress, descriptors, potential)

#         m = size(Aei,2) # total number of descriptors
#         if i == 1            
#             ce = ce*zeros(Float64,m,n)
#             cf = cf*zeros(Float64,m,n)
#             cef = cef*zeros(Float64,m,n)
#             cefs = cefs*zeros(Float64,m,n)
#             De = De*zeros(m,m,n)
#             Df = Df*zeros(m,m,n)
#             Ds = Ds*zeros(m,m,n)
#             de = de*zeros(m,n)
#             df = df*zeros(m,n)
#             ds = ds*zeros(m,n)
#             Le = Le*zeros(m,m)
#             Lf = Lf*zeros(m,m)
#             Ls = Ls*zeros(m,m)
#             he = he*zeros(m)
#             hf = hf*zeros(m)
#             hs = hs*zeros(m)
#         end
        
#         De[:,:,i] = applytranspose(Aei, Aei, m, m)    
#         Df[:,:,i] = applytranspose(Afi, Afi, m, m)    
#         Ds[:,:,i] = applytranspose(Asi, Asi, m, m)    
#         de[:,i] = applytranspose(Aei, bei, m, 1)    
#         df[:,i] = applytranspose(Afi, bfi, m, 1)    
#         ds[:,i] = applytranspose(Asi, bsi, m, 1)    

#         Aei, bei = applyweight(Aei, bei, config.we)
#         Afi, bfi = applyweight(Afi, bfi, config.wf)
#         Asi, bsi = applyweight(Asi, bsi, config.ws)
                
#         Le += applytranspose(Aei, Aei, m, m)    
#         Lf += applytranspose(Afi, Afi, m, m)    
#         Ls += applytranspose(Asi, Asi, m, m)    
#         he += applytranspose(Aei, bei, m, 1)    
#         hf += applytranspose(Afi, bfi, m, 1)    
#         hs += applytranspose(Asi, bsi, m, 1)            
        
#         if length(data)==1
#             if optim.method == "lsq"
#                 cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
#             else
#                 cei, cfi, cefi, cefsi = linearsolve(Aei, Afi, Asi, bei, bfi, bsi)                        
#             end            
#             ce = insertarray(ce, cei, i)
#             cf = insertarray(cf, cfi, i)
#             cef = insertarray(cef, cefi, i)
#             cefs = insertarray(cefs, cefsi, i)

#             emae[i,1] = calmae(Aei, bei, cei)
#             emae[i,2] = calmae(Aei, bei, cfi)
#             emae[i,3] = calmae(Aei, bei, cefi)
#             emae[i,4] = calmae(Aei, bei, cefsi)
#             fmae[i,1] = calmae(Afi, bfi, cei)
#             fmae[i,2] = calmae(Afi, bfi, cfi)
#             fmae[i,3] = calmae(Afi, bfi, cefi)
#             fmae[i,4] = calmae(Afi, bfi, cefsi)
#             smae[i,1] = calmae(Asi, bsi, cei)
#             smae[i,2] = calmae(Asi, bsi, cfi)
#             smae[i,3] = calmae(Asi, bsi, cefi)
#             smae[i,4] = calmae(Asi, bsi, cefsi)
            
#             return ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
#         end

#         if (i == 1) | (optim.method == "lsq")
#             Ae = 1.0*Aei
#             Af = 1.0*Afi
#             As = 1.0*Asi
#             be = 1.0*bei
#             bf = 1.0*bfi
#             bs = 1.0*bsi
#         else         
#             Ae, be = catarray(Ae, be, Aei, bei)
#             Af, bf = catarray(Af, bf, Afi, bfi)
#             As, bs = catarray(As, bs, Asi, bsi)                        
#         end                        

#         if (optim.method == "lsq")
#             #display(Le)
#             cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
#         else
#             cei, cfi, cefi, cefsi = linearsolve(Ae, Af, As, be, bf, bs)
#         end

#         # if i==n
#         #     display(i)
#         #     A = Le+Lf+Ls
#         #     b = he+hf+hs                
#         #     fileID = open("Ajl","w"); write(fileID,A[:]); close(fileID);
#         #     fileID = open("bjl","w"); write(fileID,b[:]); close(fileID);
#         #     display(A)    
#         #     display(b)    
#         #     coeff = (0.5*(A+A'))\b
#         #     display(coeff)
#         #     # display(size(Ae))
#         #     # display(size(Af))
#         #     # display(size(As))
#         #     coeff = cat(cat(Ae, Af, dims=1), As, dims=1)\[be[:]; bf[:]; bs[:]]
#         #     # A = Ae+Af+As
#         #     # b = be+bf+bs                        
#         #     # coeff = A\b
#         #     display(coeff)
#         #     error("here")
#         # end

#         ce = insertarray(ce, cei, i)
#         cf = insertarray(cf, cfi, i)
#         cef = insertarray(cef, cefi, i)
#         cefs = insertarray(cefs, cefsi, i)

#         emae[i,1] = calmae(Ae, be, cei)
#         emae[i,2] = calmae(Ae, be, cfi)
#         emae[i,3] = calmae(Ae, be, cefi)
#         emae[i,4] = calmae(Ae, be, cefsi)
#         fmae[i,1] = calmae(Af, bf, cei)
#         fmae[i,2] = calmae(Af, bf, cfi)
#         fmae[i,3] = calmae(Af, bf, cefi)
#         fmae[i,4] = calmae(Af, bf, cefsi)
#         smae[i,1] = calmae(As, bs, cei)
#         smae[i,2] = calmae(As, bs, cfi)
#         smae[i,3] = calmae(As, bs, cefi)
#         smae[i,4] = calmae(As, bs, cefsi)
#     end

#     if lossfunc == "energy" 
#         coeff = ce[:,end]
#     elseif lossfunc == "force" 
#         coeff = cf[:,end]
#     elseif lossfunc == "energyforce" 
#         coeff = cef[:,end]
#     elseif lossfunc == "energyforcestress" 
#         coeff = cefs[:,end]
#     end    
    
#     return coeff, ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
# end

