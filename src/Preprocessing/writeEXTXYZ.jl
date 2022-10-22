function writeEXTXYZ(filename,config,species)

    nat = [0; cumsum(config.natom[:])];
    n = config.nconfigs
    open(filename,"w") do io
        for i = 1:n
            natom = config.natom[i]            
            lattice = config.lattice[:,i]
            l1 = lattice[1] 
            l2 = lattice[2] 
            l3 = lattice[3] 
            l4 = lattice[4] 
            l5 = lattice[5] 
            l6 = lattice[6] 
            l7 = lattice[7] 
            l8 = lattice[8] 
            l9 = lattice[9] 
            energy = config.e[i]
            stress = config.stress[:,i]
            s1 = stress[1] 
            s2 = stress[2] 
            s3 = stress[3] 
            s4 = stress[4] 
            s5 = stress[5] 
            s6 = stress[6] 
            s7 = stress[7] 
            s8 = stress[8] 
            s9 = stress[9] 
            t = config.t[(nat[i]+1):nat[i+1]]
            x = config.x[:,(nat[i]+1):nat[i+1]]
            f = config.f[:,(nat[i]+1):nat[i+1]]

            line1 = """Lattice="$l1 $l2 $l3 $l4 $l5 $l6 $l7 $l8 $l9" Properties=species:S:1:pos:R:3:forces:R:3"""
            line2 = """ energy=$energy stress="$s1 $s2 $s3 $s4 $s5 $s6 $s7 $s8 $s9" pbc="T T T"\n"""
            line = line1 * line2
            println(io,natom)    
            write(io,line)
            for j = 1:natom
                x1 = x[1,j]
                x2 = x[2,j]
                x3 = x[3,j]
                f1 = f[1,j]
                f2 = f[2,j]
                f3 = f[3,j]                
                line = species[t[j]] * """  $x1  $x2  $x3  $f1  $f2  $f3\n"""
                write(io,line)
            end    
        end            
    end    
end

function convert2EXTXYZ(datapath, dataformat, fileextension, species)

    pbc = nothing 
    original = 1
    config = readconfigpath(datapath, dataformat, fileextension, species, nothing, pbc, original)
    ii = findlast("/", datapath);        
    #filename = datapath * datapath[ii[end]:end] * ".xyz"
    filename = datapath[1:ii[end]] * datapath[ii[end]:end] * ".xyz"
    writeEXTXYZ(filename, config, species)
end

function convert2EXTXYZ(data)
    for n = 1:length(data)
        datapath, dataformat, fileextension, percentage, randomize, species = getdata(data[n])
        convert2EXTXYZ(datapath, dataformat, fileextension, species)
        display(datapath)
    end
end

