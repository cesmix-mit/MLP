function readJSONpymatgen(filename, species)

    config = initializeconfig(0);  

    mystr = read(filename, String);
    data = JSON.parse(mystr);
    nconfigs = length(data)

    natoms = zeros(Int64,1,nconfigs)
    energies = zeros(1,nconfigs)
    stresses = zeros(6, nconfigs)
    forces = []
    coords = []
    types = []
    a = zeros(3, nconfigs)
    b = zeros(3, nconfigs)
    c = zeros(3, nconfigs)
    latvol = zeros(nconfigs)
    for i = 1:nconfigs
        energies[i] = data[i]["outputs"]["energy"]
        stresses[:,i] = data[i]["outputs"]["virial_stress"]
        tm = data[i]["outputs"]["forces"]
        natoms[i] = length(tm)
        f = zeros(3,natoms[i])
        for j = 1:natoms[i]
            f[:,j] = tm[j]
        end
        if i == 1
            forces = 1.0*f
        else
            forces = [forces f]
        end
        latvol[i] = data[i]["structure"]["lattice"]["volume"]
        a[:,i] = data[i]["structure"]["lattice"]["matrix"][1]
        b[:,i] = data[i]["structure"]["lattice"]["matrix"][2]
        c[:,i] = data[i]["structure"]["lattice"]["matrix"][3]
        s = data[i]["structure"]["sites"]
        t = zeros(Int64, 1, natoms[i])
        for j = 1:natoms[i]
            f[:,j] = s[j]["xyz"]
            elem = s[j]["species"][1]["element"]   
            for k = 1:length(species)
                if elem == species[k]
                    t[j] = k
                end
            end
        end
        if i == 1
            coords = 1.0*f
            types = 1*t
        else
            coords = [coords f]
            types = [types t]
        end
    end

    config.e = energies 
    config.nce = size(config.e,1);

    config.stress = stresses
    config.ncs = size(config.stress,1);
    
    config.x = coords
    config.ncx = size(config.x,1);
    config.dim = size(config.x,1);

    config.f = forces 
    config.ncf = size(config.f,1);

    config.t = types 
    config.a = a
    config.b = b
    config.c = c

    config.lattice = vcat(vcat(a, b), c)    
    config.natom = natoms
    config.nconfigs = nconfigs;
    config.natomall = sum(config.natom); 

    # display(config.lattice)
    # display(config.x)

    # dim = 3
    # if size(config.lattice,1) == dim*dim    
    #     config.a = config.lattice[1:dim,:];
    #     config.b = config.lattice[(dim+1):2*dim,:];
    #     if dim==3
    #         config.c = config.lattice[(2*dim+1):3*dim,:];
    #     end
    #     n = size(config.lattice,2)
    #     natom = [0; cumsum(config.natom[:])];

    #     # rotate lattice and positions 
    #     for i = 1:n
    #         A = config.a[:,i]
    #         B = config.b[:,i]
    #         C = config.c[:,i]
    #         a, b, c = rotatelattice(A, B, C)
    #         config.a[:,i] = a
    #         config.b[:,i] = b
    #         config.c[:,i] = c
    #         X = config.x[:,(natom[i]+1):natom[i+1]]
    #         x, R = transformcoords(a,b,c,A,B,C,X)
    #         config.x[:,(natom[i]+1):natom[i+1]] = x
    #     end
    # end
    #config.lattice = vcat(vcat(config.a, config.b), config.c)

    return config  
end






