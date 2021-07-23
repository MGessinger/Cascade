module Jade

using Nemo

export diffOp,jade_ode,acb_ode,arb_ode,monodromy,powerSeries,setPolynomial!,jade_ode_legendre,jade_ode_bessel,jade_ode_hypgeom

abstract type diffOp end

mutable struct jade_ode{P<:PolyElem} <: diffOp
	order::Int
	polys::Vector{P}
	odeC::Ptr{Nothing}
end

global const arb_ode = jade_ode{arb_poly}
global const acb_ode = jade_ode{acb_poly}

function __init__()
	print("\n")
	print("  |  _  __  ,__ \n")
	print("  | |_| | \\ |__ \n")
	print("\\_/ | | |_/ |__ \n")
	print("\nThis is Jade v.1.0, an interface to Cascade,\n\n")
	print("The C-Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations!\n")
end

# Constructors

function jade_ode{P}(polys::Vector{P}) where P <: PolyElem
	if (length(polys) == 0)
		return
	end
	deg = 0
	order = length(polys)
	while iszero(polys[order])
		order -= 1
		if order < 0
			return
		end
	end
	A = jade_ode{P}(order-1,polys[1:order],C_NULL)
	finalizer(deleteC, A)
	return A
end

function jade_ode(polys::Vector{P}) where P <: PolyElem
	return jade_ode{P}(polys)
end

# In-Place Modifications

function setPolynomial!(ode::jade_ode{T}, index::Integer, polynomial::T) where T <: PolyElem
	if ode.odeC != C_NULL
		deleteC(ode)
	end
	if index > ode.order+1
		if iszero(polynomial)
			return
		end
		arr = Vector{T}(undef,index)
		for i = 1:index-1
			if i <= ode.order+1
				arr[i] = ode.polys[i]
			else
				arr[i] = zero(polynomial)
			end
		end
		ode.polys = arr
		ode.order = index-1
	end
	ode.polys[index] = polynomial
	# The number of polynomials might have changed
	while iszero(ode.polys[end])
		ode.order -= 1
		pop!(ode.polys)
	end
	return
end

# C-interface

function translateC(ode::acb_ode)
	if ode.odeC != C_NULL
		deleteC(ode)
	end
	if ode.order == 0
		ode.odeC = C_NULL
		return C_NULL
	end
	deg = 0
	for p in ode.polys
		if deg < degree(p)
			deg = degree(p)
		end
	end
	prec = ode.polys[i].parent.base_ring.prec
	A = ccall((:acb_ode_init_blank,"libcascade"), Ptr{Nothing}, (Cint,Cint), deg, ode.order)
	for i = 1:ode.order+1
		q = Ref{acb_poly}(acb_poly(ode.polys[i],prec))
		ccall((:acb_ode_set_poly,:libcascade),Cvoid,(Ptr{Nothing},Ref{acb_poly},Cint),A,q,i-1)
	end
	ode.odeC = A
	return A
end

function deleteC(ode::diffOp)
	if ode.odeC == C_NULL
		return
	end
	ccall((:acb_ode_clear,"libcascade"), Cvoid, (Ptr{Nothing},), ode.odeC)
	ode.odeC = C_NULL
	return
end

# Cacsade interface

function powerSeries(ode::jade_ode{T},p::T,n::Integer) where T <: Union{acb_poly,arb_poly}
	if n <= degree(p)
		return p
	end
	if ode.odeC == C_NULL
		translateC(ode)
	end
	R = ode.polys[1].parent
	prec = R.base_ring.prec
	p2 = acb_poly(p,prec)
	q = Ref{acb_poly}(p2)
	ccall((:acb_ode_solve_fuchs,"libcascade"),Cvoid,(Ref{acb_poly},Ptr{Nothing},Cint,Cint),q,ode.odeC,n,prec)
	if isa(ode,acb_ode) || isa(p,acb_poly)
		return R(p2)
	else
		R2 = AcbPolyRing(ComplexField(prec),:z)
		p2 = R2(p2)
		# This is a bit of a hack, because Nemo does not support demotion of strictly real acb_polys
		return R([real(coeff(p2,i)) for i = 0:degree(p2)])
	end
end

function monodromy(ode::acb_ode,z0=0)
	if ode.odeC == C_NULL
		translateC(ode)
	end
	R = ode.polys[1].parent.base_ring
	S = MatrixSpace(R,ode.order,ode.order)
	mono = S(1)
	M = Ref{acb_mat}(mono)
	Z = Ref{acb}(R(z0))
	ccall((:find_monodromy_matrix,"libcascade"),Cvoid,(Ref{acb_mat},Ptr{Nothing},Ref{acb},Cint),M,ode.odeC,Z,R.prec)
	mono.base_ring = R
	return mono
end

function monodromy(ode::jade_ode,z0=0)
	# This function is only called iff ode is *not* an acb_ode
	if ode.odeC == C_NULL
		translateC(ode)
	end
	@warn "Monodromy matrices cannot in general be computed exactly, even for integer initial conditions!\nUsing an internal precision of 64 bits!"
	R = ComplexField(64)
	S = MatrixSpace(R,ode.order,ode.order)
	mono = S(1)
	M = Ref{acb_mat}(mono)
	Z = Ref{acb}(R(z0))
	ccall((:find_monodromy_matrix,"libcascade"),Cvoid,(Ref{acb_mat},Ptr{Nothing},Ref{acb},Cint),M,ode.odeC,Z,64)
	mono.base_ring = R
	return mono
end

# External code files come here
include("examples.jl")
include("arith.jl")

end
