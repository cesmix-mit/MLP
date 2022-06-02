function readANN(filename, species)

    config = initializeconfig(0);  

    f = readlines(filename);
    i1 = findlast("=",f[1])
    i2 = findlast("eV",f[1])
    s = f[1][(i1[1]+2):(i2[1]-2)]
    energy = parse(Float64,s) #+ 1.781565308153280e+04

    n = 0
    for j = 2:100
        if f[j] == "PRIMVEC"
            n = j+1
            break;
        end
    end
    
    a = split(f[n], ' ', keepempty=false);
    b = split(f[n+1], ' ', keepempty=false);
    c = split(f[n+2], ' ', keepempty=false);
    nat = split(f[n+4], ' ', keepempty=false);

    natom = parse(Int32,nat[1])
    a = parse.(Float64,a)
    b = parse.(Float64,b)
    c = parse.(Float64,c)

    m = n+4
    s = String.(reduce(vcat, split.(f[(m+1):(m+natom)])))
    s = reshape(s, (Int64(length(s)/natom), natom))
    
    frame = parse.(Float64, s[2:end,:]); 

    t = ones(Int32,1,natom);
    for k = 1:length(species)
        ind = findall(s[1,:] .== species[k]);
        t[ind] .= k;
    end              

    config.x = frame[1:3,:]
    config.f = frame[4:6,:]
    config.e = reshape([energy],(1,1));
    config.t = t
    config.lattice = reshape([a[:]; b[:]; c[:]], (9,1))

    config.a = reshape(a,(3,1))
    config.b = reshape(b,(3,1))
    config.c = reshape(c,(3,1))

    config.nce = size(config.e,1);    
    config.ncx = size(config.x,1);
    config.dim = size(config.x,1);
    config.ncf = size(config.f,1);
    config.natom = reshape([natom],(1,1))
    config.nconfigs = 1;
    config.natomall = sum(config.natom); 

    return config 
end

