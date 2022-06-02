function istricliniclattice(a1, a2, a3)

if (abs(a1[2])<1e-12) & (abs(a1[3])<1e-12) & (abs(a2[3])<1e-12)
    return true
else
    return false
end

end

function tricliniclatticevectors(A1, A2, A3)

A1norm = norm(A1);
A2norm = norm(A2);
A3norm = norm(A3);
A1hat = A1/A1norm;

lx = A1norm;
xy = dot(A2,A1hat);
ly = sqrt(A2norm^2 - xy^2); #norm(cross(A1hat,A2));
xz = dot(A3,A1hat);
yz = (dot(A2, A3) - xy*xz)/ly;
lz = sqrt(A3norm^2 - xz^2 - yz^2);

a1 = [lx; 0; 0];
a2 = [xy; ly; 0];
a3 = [xz; yz; lz];

return a1, a2, a3, lx, ly, lz, xy, xz, yz  

end

function tricliniclatticevectors(a, b, c, cosalpha, cosbeta, cosgamma)

lx = a
xy = b*cosgamma;
xz = c*cosbeta;
ly = sqrt(b*b - xy*xy);
yz = (b*c*cosalpha - xy*xz)/ly;
lz = sqrt(c*c - xz*xz - yz*yz);

a1 = [lx; 0; 0];
a2 = [xy; ly; 0];
a3 = [xz; yz; lz];

return a1, a2, a3, lx, ly, lz, xy, xz, yz  

end

function triclinicbox(a, b, c, cosalpha, cosbeta, cosgamma)

    lx = 1.0*a
    xy = b.*cosgamma;
    xz = c.*cosbeta;
    ly = sqrt.(b.*b .- xy.*xy);
    yz = (b.*c.*cosalpha .- xy.*xz)./ly;
    lz = sqrt.(c.*c .- xz.*xz .- yz.*yz);
        
    if (length(lz) > length(lx)) 
        lx = a*ones(size(lz))
    end

    return lx, ly, lz, xy, xz, yz     
end


function tricliniclatticeparameters(lx, ly, lz, xy, xz, yz)

a = lx
b = sqrt(ly*ly + xy*xy)
c = sqrt(lz*lz + xz*xz + yz*yz)
cosalpha = (xy*xz + ly*yz)/(b*c);
cosbeta = xz/c;
cosgamma = xy/b;

return a, b, c, cosalpha, cosbeta, cosgamma

end

function tricliniclatticeparameters(a1, a2, a3)

a, b, c, cosalpha, cosbeta, cosgamma = tricliniclatticeparameters(a1[1], a2[2], a3[3], a2[1], a3[1], a3[2])

return a, b, c, cosalpha, cosbeta, cosgamma

end

function tricliniclattice(A1, A2, A3)

a1, a2, a3 = tricliniclatticevectors(A1, A2, A3)
a, b, c, cosalpha, cosbeta, cosgamma = triclinicparameters(a1[1], a2[2], a3[3], a2[1], a3[1], a3[2])

return a, b, c, cosalpha, cosbeta, cosgamma, a1, a2, a3

end
    
function tricliniclattice(lx, ly, lz, xy, xz, yz, fraccoords)

    N = length(lx)
    dim, M = size(fraccoords)
    AtomBasis = zeros(dim,M,N)
    for n = 1:N
        for m = 1:M
            AtomBasis[1,m,n] = lx[n]*fraccoords[1,m] + xy[n]*fraccoords[2,m] + xz[n]*fraccoords[3,m]
            AtomBasis[2,m,n] = ly[n]*fraccoords[2,m] + yz[n]*fraccoords[3,m]      
            AtomBasis[3,m,n] = lz[n]*fraccoords[3,m]             
            # AtomBasis[1,m,n] = lx[n]*fraccoords[1,m]
            # AtomBasis[2,m,n] = xy[n]*fraccoords[1,m] + ly[n]*fraccoords[2,m]       
            # AtomBasis[3,m,n] = xz[n]*fraccoords[1,m] + yz[n]*fraccoords[2,m] + lz[n]*fraccoords[3,m]             
        end
    end

    A1 = zeros(3, N)
    A2 = zeros(3, N)
    A3 = zeros(3, N)
    A1[1,:] = lx
    A2[1,:] = xy
    A2[2,:] = ly
    A3[1,:] = xz
    A3[2,:] = yz    
    A3[3,:] = lz

    return A1, A2, A3, AtomBasis
end








