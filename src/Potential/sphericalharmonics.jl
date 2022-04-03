# Copy from Christoph Ortner's ACE.jl/src/polynomials/sphericalharmonics.jl

using Base: @kwdef

@kwdef mutable struct SphericalCoords
	r::Float64 = 1.0
	cosφ::Float64 = 1.0
	sinφ::Float64 = 0.0
	cosθ::Float64 = 1.0
	sinθ::Float64 = 0.0
end

spher2cart(S::SphericalCoords) = S.r * [S.cosφ*S.sinθ, S.sinφ*S.sinθ, S.cosθ]

function cart2spher(R::Vector{Float64})
	@assert length(R) == 3
	r = norm(R)
	φ = atan(R[2], R[1])
	sinφ, cosφ = sincos(φ)
	cosθ = R[3] / r
	sinθ = sqrt(R[1]^2+R[2]^2) / r
	return SphericalCoords(r, cosφ, sinφ, cosθ, sinθ)
end

function cart2spher!(S::SphericalCoords, R::Vector{Float64})
	@assert length(R) == 3
	S.r = sqrt(R[1]*R[1] + R[2]*R[2]+ R[3]*R[3])
	φ = atan(R[2], R[1])
	S.sinφ, S.cosφ = sincos(φ)
	S.cosθ = R[3] / S.r
	S.sinθ = sqrt(R[1]^2+R[2]^2) / S.r
end

SphericalCoords(φ, θ) = SphericalCoords(1.0, cos(φ), sin(φ), cos(θ), sin(θ))
SphericalCoords(r, φ, θ) = SphericalCoords(r, cos(φ), sin(φ), cos(θ), sin(θ))

function dspher_to_dcart(S, f_φ_div_sinθ, f_θ)
	r = S.r + eps()

    g = [- (S.sinφ * f_φ_div_sinθ) + (S.cosφ * S.cosθ * f_θ), 
            (S.cosφ * f_φ_div_sinθ) + (S.sinφ * S.cosθ * f_θ), 
            - ( S.sinθ * f_θ) ] / r

   return g
end

sizeP(maxL) = div((maxL + 1) * (maxL + 2), 2)
sizeY(maxL) = (maxL + 1) * (maxL + 1)
indexP(l::Integer,m::Integer) = m + div(l*(l+1), 2) + 1
indexY(l::Integer, m::Integer) = m + l + (l*l) + 1

function idx2lm(i::Integer) 
	l = floor(Int, sqrt(i-1) + 1e-10)
	m = i - (l + (l*l) + 1)
	return l, m 
end 

struct ALCoefficients
	L::Int
	A::Vector{Float64}
	B::Vector{Float64}
end

function ALCoefficients(L::Integer)
	# Precompute coefficients ``a_l^m`` and ``b_l^m`` for all l <= L, m <= l
	alp = ALCoefficients(L, zeros(sizeP(L)), zeros(sizeP(L)))
	for l in 2:L
		ls = l*l
		lm1s = (l-1) * (l-1)
		for m in 0:(l-2)
			ms = m * m
			alp.A[indexP(l, m)] = sqrt((4 * ls - 1.0) / (ls - ms))
			alp.B[indexP(l, m)] = -sqrt((lm1s - ms) / (4 * lm1s - 1.0))
		end
	end
	return alp
end

function ALPolynomials(P, alp::ALCoefficients, S::SphericalCoords)
	L = alp.L 
	A = alp.A 
	B = alp.B 
	@assert length(A) >= sizeP(L)
	@assert length(B) >= sizeP(L)
	@assert length(P) >= sizeP(L)

	temp = sqrt(0.5/π)
	P[indexP(0, 0)] = temp
	if L == 0; return P; end

	P[indexP(1, 0)] = S.cosθ * sqrt(3) * temp
	temp = - sqrt(1.5) * S.sinθ * temp
	P[indexP(1, 1)] = temp

	for l in 2:L
		il = ((l*(l+1)) ÷ 2) + 1
		ilm1 = il - l
		ilm2 = ilm1 - l + 1
		for m in 0:(l-2)
			@inbounds P[il+m] = A[il+m] * (     S.cosθ * P[ilm1+m]
  					                           + B[il+m] * P[ilm2+m] )
		end
		@inbounds P[il+l-1] = S.cosθ * sqrt(2 * (l - 1) + 3) * temp
		temp = -sqrt(1.0 + 0.5 / l) * S.sinθ * temp
		@inbounds P[il+l] = temp
	end

	return P
end

function ALPolynomials(P, dP, alp::ALCoefficients, S::SphericalCoords)
	L = alp.L 
	A = alp.A 
	B = alp.B 
	@assert length(A) >= sizeP(L)
	@assert length(B) >= sizeP(L)
	@assert length(P) >= sizeP(L)
	@assert length(dP) >= sizeP(L)

	temp = sqrt(0.5/π)
	P[indexP(0, 0)] = temp
	temp_d = 0.0
	dP[indexP(0, 0)] = temp_d
	if L == 0; return P, dP; end

	P[indexP(1, 0)] = S.cosθ * sqrt(3) * temp
	dP[indexP(1, 0)] = -S.sinθ * sqrt(3) * temp + S.cosθ * sqrt(3) * temp_d
	temp1, temp_d = ( - sqrt(1.5) * temp,
					      - sqrt(1.5) * (S.cosθ * temp + S.sinθ * temp_d) )
	P[indexP(1, 1)] = temp1
	dP[indexP(1, 1)] = temp_d

	for l in 2:L
		m = 0
		@inbounds P[indexP(l, m)] =
				A[indexP(l, m)] * (     S.cosθ * P[indexP(l - 1, m)]
				             + B[indexP(l, m)] * P[indexP(l - 2, m)] )
		@inbounds dP[indexP(l, m)] =
			A[indexP(l, m)] * (
							- S.sinθ * P[indexP(l - 1, m)]
							+ S.cosθ * dP[indexP(l - 1, m)]
			             + B[indexP(l, m)] * dP[indexP(l - 2, m)] )

		for m in 1:(l-2)
			@inbounds P[indexP(l, m)] =
					A[indexP(l, m)] * (     S.cosθ * P[indexP(l - 1, m)]
					             + B[indexP(l, m)] * P[indexP(l - 2, m)] )
			@inbounds dP[indexP(l, m)] =
				A[indexP(l, m)] * (
								- S.sinθ^2 * P[indexP(l - 1, m)]
								+ S.cosθ * dP[indexP(l - 1, m)]
				             + B[indexP(l, m)] * dP[indexP(l - 2, m)] )
		end
		@inbounds P[indexP(l, l - 1)] = sqrt(2 * (l - 1) + 3) * S.cosθ * temp1
		@inbounds dP[indexP(l, l - 1)] = sqrt(2 * (l - 1) + 3) * (
									        -S.sinθ^2 * temp1 + S.cosθ * temp_d )

      (temp1, temp_d) = (
					-sqrt(1.0+0.5/l) * S.sinθ * temp1,
		         -sqrt(1.0+0.5/l) * (S.cosθ * temp1 * S.sinθ + S.sinθ * temp_d) )
		@inbounds P[indexP(l, l)] = temp1
		@inbounds dP[indexP(l, l)] = temp_d
	end

	return P, dP
end

struct SHBasis
    L::Int
	alp::ALCoefficients
	S::SphericalCoords
	P::Vector{Float64}
    dP::Vector{Float64}    
    Y::Vector{Complex{Float64}}
    dY::Matrix{Complex{Float64}}
    cg::Vector{Float64}    
    indl::Matrix{Int32}     
    indm::Matrix{Int32}         
    rowm::Vector{Int32}        
    #dY::Vector{Vector{Complex{Float64}}}
end

function SHBasis(L::Integer)    
    alp = ALCoefficients(L)
    S = SphericalCoords(1.0, 0.0, 0.0)
    P = zeros(sizeP(L))
    dP = zeros(sizeP(L))
    Y = zeros(Complex{Float64}, sizeY(L))
    dY = zeros(Complex{Float64}, sizeY(L), 3)
    indl = cgconditions(L);
    cg,indm,rowm = cgcoefficients(indl);        
    return SHBasis(L, alp, S, P, dP, Y, dY, cg, indl, indm, rowm)
end

function cYlm!(Y, L, S::SphericalCoords, P)
	@assert length(P) >= sizeP(L)
	@assert length(Y) >= sizeY(L)
    @assert abs(S.cosθ) <= 1.0

	ep = 1 / sqrt(2) + im * 0
	for l = 0:L
		Y[indexY(l, 0)] = P[indexP(l, 0)] * ep
	end

   sig = 1
   ep_fact = S.cosφ + im * S.sinφ
	for m in 1:L
		sig *= -1
		ep *= ep_fact            # ep =   exp(i *   m  * φ)
		em = sig * conj(ep)      # ep = ± exp(i * (-m) * φ)
		for l in m:L
			p = P[indexP(l,m)]
			@inbounds Y[indexY(l, -m)] = em * p   # (-1)^m * p * exp(-im*m*phi) / sqrt(2)
			@inbounds Y[indexY(l,  m)] = ep * p   #          p * exp( im*m*phi) / sqrt(2)
		end
	end

	return Y
end

function cYlm!(Y, L, S::SphericalCoords, P, alp::ALCoefficients, R::Vector{Float64})
	cart2spher!(S, R)
	ALPolynomials(P, alp, S)
	cYlm!(Y, L, S, P)
	return Y 
end

function cYlm!(SH::SHBasis, R::Vector{Float64})
	cart2spher!(SH.S, R)
	ALPolynomials(SH.P, SH.alp, SH.S)
	cYlm!(SH.Y, SH.L, SH.S, SH.P)
	return SH 
end

function cYlm!(SH::SHBasis, R::Matrix{Float64})
    dim, N = size(R)
    Ylm = zeros(Complex{Float64}, N, sizeY(SH.L))
    for n = 1:N
        cart2spher!(SH.S, R[:,n])
        ALPolynomials(SH.P, SH.alp, SH.S)
        cYlm!(SH.Y, SH.L, SH.S, SH.P)
        Ylm[n,:] = SH.Y 
    end
	return Ylm 
end

function cYlm_de!(Y, dY, L, S::SphericalCoords, P, dP)
	@assert length(P) >= sizeP(L)
	@assert length(Y) >= sizeY(L)
	@assert size(dY,1) >= sizeY(L)

	# m = 0 case
	ep = 1 / sqrt(2)
	for l = 0:L
		Y[indexY(l, 0)] = P[indexP(l, 0)] * ep 
		#dY[indexY(l, 0)] = dspher_to_dcart(S, 0.0, dP[indexP(l, 0)] * ep)
        dY[indexY(l, 0),:] = dspher_to_dcart(S, 0.0, dP[indexP(l, 0)] * ep) 
	end

   sig = 1
   ep_fact = S.cosφ + im * S.sinφ

	for m in 1:L
		sig *= -1
		ep *= ep_fact            # ep =   exp(i *   m  * φ)
		em = sig * conj(ep)      # ep = ± exp(i * (-m) * φ)
		dep_dφ = im *   m  * ep
		dem_dφ = im * (-m) * em
		for l in m:L
			p_div_sinθ = P[indexP(l,m)]
			@inbounds Y[indexY(l, -m)] = em * p_div_sinθ * S.sinθ 
			@inbounds Y[indexY(l,  m)] = ep * p_div_sinθ * S.sinθ 

			dp_dθ = dP[indexP(l,m)]
			@inbounds dY[indexY(l, -m),:] = dspher_to_dcart(S, dem_dφ * p_div_sinθ, em * dp_dθ) 
			@inbounds dY[indexY(l,  m),:] = dspher_to_dcart(S, dep_dφ * p_div_sinθ, ep * dp_dθ) 
		end
	end

	return Y, dY
end

function cYlm_de!(Y, dY, L, S::SphericalCoords, P, dP, alp::ALCoefficients, R::Vector{Float64})
	cart2spher!(S, R)
	ALPolynomials(P, dP, alp, S)
	cYlm_de!(Y, dY, L, S, P, dP)
	return Y, dY 
end

function cYlm_de!(SH::SHBasis, R::Vector{Float64})
	cart2spher!(SH.S, R)
	ALPolynomials(SH.P, SH.dP, SH.alp, SH.S)
	cYlm_de!(SH.Y, SH.dY, SH.L, SH.S, SH.P, SH.dP)
	return SH 
end

function cYlm_de!(SH::SHBasis, R::Matrix{Float64})
    dim, N = size(R)
    M = sizeY(SH.L)
    Ylm = zeros(Complex{Float64}, N, M)
    dYlm = zeros(Complex{Float64}, N, M, 3)
    for n = 1:N        
        cart2spher!(SH.S, R[:,n])
        ALPolynomials(SH.P, SH.dP, SH.alp, SH.S)
        cYlm_de!(SH.Y, SH.dY, SH.L, SH.S, SH.P, SH.dP)
        Ylm[n,:] = SH.Y 
        dYlm[n,:,:] = SH.dY 
    end    
    dYlm = reshape(permutedims(dYlm, (3,1,2)), (3, N, M))
	return Ylm, dYlm    
end


# r = 1.0
# φ = pi/4.5 
# θ = pi/2.3
# lmax = 4 
# alp = ALCoefficients(lmax)
# s = SphericalCoords(r, φ, θ)
# P = zeros(sizeP(lmax))
# P = ALPolynomials(P, alp, s)
# Y = zeros(Complex{Float64}, sizeY(lmax))
# Y = cYlm!(Y, lmax, s, P)

# P = zeros(sizeP(lmax))
# dP = zeros(sizeP(lmax))
# P,dP = ALPolynomials(P, dP, alp, s)

# Y = zeros(Complex{Float64}, sizeY(lmax))
# #dY = [zeros(Complex{Float64}, 3) for _ in 1:sizeY(lmax)]
# dY = zeros(Complex{Float64}, sizeY(lmax), 3)
# Y, dY = cYlm_de!(Y, dY, lmax, s, P, dP)

# using SphericalHarmonics
# P = zeros(sizeP(lmax))
# P = ALPolynomials(P, alp, s)
# P1 = computePlmcostheta(θ, lmax)
# Y1 = computeYlm(θ, φ, lmax)

# display(P - P1)
# display(Y - Y1)

# R = spher2cart(s)
# sh = SHBasis(lmax)
# sh = cYlm!(sh, R)
# sh = cYlm_de!(sh, R)

# rij = rand(3,10)
# Ylm = cYlm!(sh, rij)
# Ylm2, dYlm2 = cYlm_de!(sh, rij)



