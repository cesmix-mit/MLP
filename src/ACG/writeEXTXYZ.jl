function writeEXTXYZ2(filename, species, natoms, lattices, atomtypes, atompos)

    nat = [0; cumsum(natoms[:])];
    n = length(natoms)
    open(filename,"w") do io
        for i = 1:n
            natom = natoms[i]            
            lattice = lattices[:,i]
            l1 = lattice[1] 
            l2 = lattice[2] 
            l3 = lattice[3] 
            l4 = lattice[4] 
            l5 = lattice[5] 
            l6 = lattice[6] 
            l7 = lattice[7] 
            l8 = lattice[8] 
            l9 = lattice[9] 
            t = atomtypes[(nat[i]+1):nat[i+1]]
            x = atompos[:,(nat[i]+1):nat[i+1]]

            line1 = """Lattice="$l1 $l2 $l3 $l4 $l5 $l6 $l7 $l8 $l9" Properties=species:S:1:pos:R:3"""
            line2 = """ pbc="T T T"\n"""
            line = line1 * line2
            println(io,natom)    
            write(io,line)
            for j = 1:natom
                x1 = x[1,j]
                x2 = x[2,j]
                x3 = x[3,j]
                line = species[t[j]] * """  $x1  $x2  $x3\n"""
                write(io,line)
            end    
        end            
    end    
end
