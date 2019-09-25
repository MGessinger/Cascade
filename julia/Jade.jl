module Jade

using Nemo

export diffOp,jade_ode,acb_ode,arb_ode,monodromy,powerSeries,setPolynomial!,jade_ode_legendre,acb_ode_bessel,acb_ode_hypgeom

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

function translateC(ode::diffOp)
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
    A = ccall((:acb_ode_init_blank,"libcascade"), Ptr{Nothing}, (Cint,Cint), deg, ode.order)
    for i = 1:ode.order+1
        q = Ref{acb_poly}(acb_poly(ode.polys[i],0))
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

function powerSeries(ode::acb_ode,p::acb_poly,n::Integer)
    if n <= degree(p)
        return p
    end
    if ode.odeC == C_NULL
        translateC(ode);
    end
    polyRing = p.parent
    q = Ref{acb_poly}(p)
    ccall((:find_power_series,"libcascade"),Cvoid,(Ref{acb_poly},Ptr{Nothing},Cint,Cint),q,ode.odeC,n,polyRing.base_ring.prec)
    return polyRing(p)
end

function monodromy(ode::diffOp,z0=0)
    if ode.odeC == C_NULL
        translateC(ode);
    end
    if isa(ode,acb_ode) || isa(ode,arb_ode)
        R = ode.polys[1].parent.base_ring
    else
        print("Plese enter a precision value: ")
        R = ComplexField(parse(Int,readline()))
    end
    S = MatrixSpace(R,ode.order,ode.order)
    mono = S(1)
    M = Ref{acb_mat}(mono)
    Z = Ref{acb}(R(z0))
    ccall((:find_monodromy_matrix,"libcascade"),Cvoid,(Ref{acb_mat},Ptr{Nothing},Ref{acb},Cint),M,ode.odeC,Z,R.prec)
    mono.base_ring = R
    return mono
end

function jade_ode_legendre(R::PolyRing,n::Integer)
    ode = jade_ode([R(n*(n+1)),R([0,-2]),R([1,0,-1])])
    ode.odeC = ccall((:acb_ode_legendre,:libcascade),Ptr{Nothing},(Cint,),n)
    return ode
end

function acb_ode_bessel(R::AcbPolyRing,nu::acb)
    p = R([0,0,1])
    ode = acb_ode([p-nu*nu,R([0,1]),R([0,0,1])])
    N = Ref{acb}(nu)
    ode.odeC = ccall((:acb_ode_bessel,:libcascade),Ptr{Nothing},(Ref{acb},Cint),N,R.base_ring.prec)
    return ode
end

function acb_ode_hypgeom(R::AcbPolyRing,a::acb,b::acb,c::acb)
    ode = acb_ode([R(-a*b),R([c,(a+b+1)]),R([0,1,-1])])
    A = Ref{acb}(a)
    B = Ref{acb}(b)
    C = Ref{acb}(c)
    ode.odeC = ccall((:acb_ode_hypgeom,:libcascade),Ptr{Nothing},(Ref{acb},Ref{acb},Ref{acb},Cint),A,B,C,R.base_ring.prec)
    return ode
end

# External code files come here
include("arith.jl")

end
