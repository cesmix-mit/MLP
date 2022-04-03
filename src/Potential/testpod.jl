using Revise, LinearAlgebra
using Base: @kwdef

include("PODdescriptors.jl")
include("bispectrum.jl")
include("radialbasis.jl")

rin = 1.0
rcut = 5.0
K = 8
L = 4
N = 100
R = rand(3,N)
rij = (rcut-rin)*R .+ rin

pod = POD(nbody=2, pdegree=[6,10], nbasis = [K], rin = rin, rcut=rcut, gamma0 = [0.0, 2, 4]);
pod.sh = SHBasis(L)
sh = SHBasis(L)

phi = radialsphericalharmonicbasis(pod, sh, rij);
A = sum(phi, dims=1);
#ps = powerspectrum(A);
ps = powerspectrum2(real(A),imag(A));
bs = bispectrum3(real(A), imag(A), sh.cg, sh.indl, sh.indm, sh.rowm);

phi2, dphi = radialsphericalharmonicbasis_de(pod, sh, rij)
display(maximum(abs.(phi[:]-phi2[:])))

phi3, dphi3 = radialsphericalharmonicbasis(pod, rij);
display(maximum(abs.(phi3[:]-phi2[:])))
display(maximum(abs.(dphi3[:]-dphi[:])))

a = pi/3; b = pi/4; c = pi/6
Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]
Rz = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1]
Rr = Rz*Ry*Rx

Rij = Rr*rij;
phi = radialsphericalharmonicbasis(pod, sh, Rij);
B = sum(phi, dims=1);
psR = powerspectrum2(real(B),imag(B));
bsR = bispectrum3(real(B), imag(B), sh.cg, sh.indl, sh.indm, sh.rowm);

display([maximum(abs.(ps[:]-psR[:])) maximum(abs.(bs[:]-bsR[:]))])

phi, dphi = radialsphericalharmonicbasis_de(pod, sh, rij)
M = size(phi,2)
dphi1 = 0.0*dphi 
for n = 1:N    
    for d = 1:3        
        rij1 = reshape(1.0*rij[:,n], (3,1)) 
        rij1[d] = rij1[d] + 1e-6
        phi1 = radialsphericalharmonicbasis(pod, sh, rij1);                      
        dphi1[d,n,:,:] = (reshape(phi1,(M,K)) - phi[n,:,:])*1e6
    end
end
display(maximum(abs.(dphi[:]-dphi1[:])))

M = size(dphi,3)
A = reshape(sum(phi, dims=1), (M,K));
ps, dps = powerspectrum1(real(A), imag(A), real(dphi), imag(dphi))
dps1 = 0.0*dps 
for n = 1:N    
    for d = 1:3        
        rij1 = 1.0*rij
        rij1[d,n] = rij1[d,n] + 1e-6
        phi1 = radialsphericalharmonicbasis(pod, sh, rij1);      
        A1 = reshape(sum(phi1, dims=1), (M,K));         
        local ps1 = powerspectrum1(real(A1), imag(A1), nothing, nothing)       
        dps1[d,n,:,:] = (ps1 - ps)*1e6
    end
end
display(maximum(abs.(dps[:]-dps1[:])))

ps, dps = powerspectrum2(real(A), imag(A), real(dphi), imag(dphi))
dps1 = 0.0*dps 
for n = 1:N    
    for d = 1:3        
        rij1 = 1.0*rij
        rij1[d,n] = rij1[d,n] + 1e-6
        phi1 = radialsphericalharmonicbasis(pod, sh, rij1);      
        A1 = reshape(sum(phi1, dims=1), (M,K));         
        local ps1 = powerspectrum2(real(A1), imag(A1), nothing, nothing)       
        dps1[d,n,:,:] = (ps1 - ps)*1e6
    end
end
display(maximum(abs.(dps[:]-dps1[:])))

#bs = bispectrum3(real(A), imag(A), cg, indl, indm, rowm);
bs, dbs = bispectrum1(real(A), imag(A), real(dphi), imag(dphi), sh.cg, sh.indl, sh.indm, sh.rowm)
dbs1 = 0.0*dbs 
for n = 1:N    
    for d = 1:3        
        rij1 = 1.0*rij
        rij1[d,n] = rij1[d,n] + 1e-6
        phi1 = radialsphericalharmonicbasis(pod, sh, rij1);      
        A1 = reshape(sum(phi1, dims=1), (M,K));         
        local bs1 = bispectrum1(real(A1), imag(A1), nothing, nothing, sh.cg, sh.indl, sh.indm, sh.rowm)       
        dbs1[d,n,:,:] = (bs1 - bs)*1e6
    end
end
display(maximum(abs.(dbs[:]-dbs1[:])))

ps, dps = powerspectrum1(real(A), imag(A), real(dphi), imag(dphi))
Ni = 1
ai = ones(Int32, N)
aj = ones(Int32, N)
Ar, Ai = radialsphericalharmonicbasis_sum(real(phi), imag(phi), ai, Ni);
ps2, dps2 = powerspectrum1(Ar, Ai, real(dphi), imag(dphi), ai, aj);
display(maximum(abs.(ps[:]-ps2[:])))
display(maximum(abs.(dps[:]-dps2[:])))

bs2, dbs2 = bispectrum1(Ar, Ai, real(dphi), imag(dphi), sh.cg, sh.indl, sh.indm, sh.rowm, ai, aj);
display(maximum(abs.(bs[:]-bs2[:])))
display(maximum(abs.(dbs[:]-dbs2[:])))

A = radialsphericalharmonicbasis_sum(phi, ai, Ni);
display(maximum(abs.(Ar[:]-real(A[:]))))
display(maximum(abs.(Ai[:]-imag(A[:]))))



