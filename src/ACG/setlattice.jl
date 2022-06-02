mutable struct LatticeStruct
    a1::Array{Float64,1};
    a2::Array{Float64,1};
    a3::Array{Float64,1};
    a::Float64
    b::Float64
    c::Float64
    cosalpha::Float64
    cosbeta::Float64
    cosgamma::Float64
    atombasis::Array{Float64,2};
    fraccoords::Array{Float64,2};
    atomtype::Array{Int64,1};    
    style::String;         
    triclinic::Bool
    LatticeStruct() = new();
end

function setlattice(style, a1=[1.0, 0.0, 0.0], a2=[0.0, 1.0, 0.0], a3=[0.0, 0.0, 1.0], atombasis=nothing, atomtype=nothing)

lat = LatticeStruct();

if (lowercase(style) == "custom") 
    if (length(a1) != 3) & (length(a2) != 3) & (length(a3) != 3)
        error("Lattice vectors must be length 3.");
    end
    if size(atombasis,1) != 3
        atombasis = atombasis';
    end
    if size(atombasis,1) != 3
        error("Atom basis is incorrect.");
    end
else
    a1 = [1.0, 0.0, 0.0];
    a2 = [0.0, 1.0, 0.0];
    a3 = [0.0, 0.0, 1.0];
    
    if (lowercase(style) === "hex") 
        a2[2] = sqrt(3.0);
    end
    if (lowercase(style) === "hcp") 
        a2[2] = sqrt(3.0);
        a3[3] = sqrt(8.0/3.0);
    end
    
    if (lowercase(style) === "sc") 
        atombasis = [0 0 0];        
    elseif (lowercase(style) === "bcc") 
        atombasis = [0 0 0; 0.5 0.5 0.5];
    elseif (lowercase(style) === "fcc") 
        atombasis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5];
    elseif (lowercase(style) === "hcp") 
        atombasis = [0 0 0; 0.5 0.5 0.0; 0.5 5.0/6.0 0.5; 0.0 1.0/3.0 0.5];
    elseif (lowercase(style) === "diamond") 
        atombasis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5; 0.25 0.25 0.25; 0.25 0.75 0.75; 0.75 0.25 0.75; 0.75 0.75 0.25];        
    elseif (lowercase(style) === "sq") 
        atombasis = [0 0 0];        
    elseif (lowercase(style) === "sq2") 
        atombasis = [0 0 0; 0.5 0.5 0.0];
    elseif (lowercase(style) === "hex") 
        atombasis = [0 0 0; 0.5 0.5 0.0];    
    end
    atombasis = atombasis'
end

lat.fraccoords = inv([a1 a2 a3])*atombasis    

triclinic = istricliniclattice(a1, a2, a3);
if triclinic == false
    display("Since the lattice is not triclinic, we convert it to triclinic system.")
    lat.a, lat.b, lat.c, lat.cosalpha, lat.cosbeta, lat.cosgamma, 
        lat.a1, lat.a2, lat.a3 = tricliniclattice(a1, a2, a3)
    lat.atombasis = [lat.a1 lat.a2 lat.a3]*lat.fraccoords    
else
    lat.a1 = a1
    lat.a2 = a2
    lat.a3 = a3
    lat.atombasis = atombasis    
    lat.a, lat.b, lat.c, lat.cosalpha, lat.cosbeta, lat.cosgamma = tricliniclatticeparameters(a1, a2, a3)
end

natombasis = size(lat.atombasis,2);
if atomtype===nothing
    lat.atomtype = Int64.(ones(natombasis));
else
    lat.atomtype = atomtype;
end

lat.style = style;
lat.triclinic = triclinic

return lat

end

function setlattice(a, b, c, alpha, beta, gamma, atombasis, atomtype=nothing)

if size(atombasis,1) != 3
    atombasis = atombasis';
end
if size(atombasis,1) != 3
    error("Atom basis is incorrect.");
end

alpha = alpha*pi/180;
beta = beta*pi/180;
gamma = gamma*pi/180;

lat = LatticeStruct();
lat.style = "custom";
lat.triclinic = true;
lat.a = a;
lat.b = b;
lat.c = c;
lat.cosalpha = cos(alpha);
lat.cosbeta = cos(beta);
lat.cosgamma = cos(gamma);
lat.atombasis = atombasis

lat.a1, lat.a2, lat.a3 = tricliniclatticevectors(lat.a, lat.b, lat.c, lat.cosalpha, lat.cosbeta, lat.cosgamma)
lat.fraccoords = inv([lat.a1 lat.a2 lat.a3])*atombasis    

natombasis = size(lat.atombasis,2);
if atomtype===nothing
    lat.atomtype = Int64.(ones(natombasis));
else
    lat.atomtype = atomtype;
end

end

function setlatticefraction(a, b, c, alpha, beta, gamma, fraccoords, atomtype=nothing)

if size(fraccoords,1) != 3
    fraccoords = fraccoords';
end
if size(fraccoords,1) != 3
    error("Atom basis is incorrect.");
end
    
lat = LatticeStruct();

alpha = alpha*pi/180;
beta = beta*pi/180;
gamma = gamma*pi/180;

cosalpha = cos(alpha);
cosbeta = cos(beta);
cosgamma = cos(gamma);
lat.a1, lat.a2, lat.a3 = tricliniclatticevectors(a, b, c, cosalpha, cosbeta, cosgamma)
atombasis = [lat.a1 lat.a2 lat.a3]*fraccoords

lat.style = "custom";
lat.triclinic = true;
lat.a = a;
lat.b = b;
lat.c = c;
lat.cosalpha = cosalpha;
lat.cosbeta = cosbeta;
lat.cosgamma = cosgamma;
lat.atombasis = atombasis
lat.fraccoords  = fraccoords 

natombasis = size(lat.atombasis,2);
if atomtype===nothing
    lat.atomtype = Int64.(ones(natombasis));
else
    lat.atomtype = atomtype;
end

return lat 

end

function referenceorigins(m, n, p)

xref = zeros(3, m*n*p)
for i = 1:p # z
    for j = 1:n # y
        for k = 1:m # x
            t = k + m*(j-1) + m*n*(i-1)
            xref[1, t] = k-1;  
            xref[2, t] = j-1;              
            xref[3, t] = i-1;              
        end
    end
end

return xref

end


function replicatefraccoords(fraccoords, m, n, p)

    if size(fraccoords,1) != 3
        fraccoords = fraccoords';
    end
    if size(fraccoords,1) != 3
        error("Atom basis is incorrect.");
    end

    xref = referenceorigins(m, n, p);

    d = 3
    nx = size(fraccoords,2)
    nl = size(xref,2)  # number of lattices

    y = zeros(d, nx, nl)
    for i = 1:nl
        for j = 1:nx                        
            y[1,j,i] = (xref[1,i] + fraccoords[1,j])/m;
            y[2,j,i] = (xref[2,i] + fraccoords[2,j])/n;
            y[3,j,i] = (xref[3,i] + fraccoords[3,j])/p;            
        end
    end
    y = reshape(y, (d, nx*nl))

    return y
end



# mutable struct LatticeArray
#     a1::Array{Float64,2};
#     a2::Array{Float64,2};
#     a3::Array{Float64,2};
#     a::Array{Float64,1}
#     b::Array{Float64,1}
#     c::Array{Float64,1}
#     cosalpha::Array{Float64,1}
#     cosbeta::Array{Float64,1}
#     cosgamma::Array{Float64,1}
#     atombasis::Array{Float64,3};
#     fraccoords::Array{Float64,3};
#     atomtype::Array{Int64,2};    
#     style::String;         
#     triclinic::boolean 
#     LatticeStruct() = new();
# end
