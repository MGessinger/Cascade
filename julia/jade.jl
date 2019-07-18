module Jade

using Nemo

export acb_ode,monodromy,powerSeries,setPolynomial,setInitialValues,acb_ode_legendre,acb_ode_bessel

mutable struct acb_ode{T<:Integer}
    order::T
    polys::Array{acb_poly,1}
    odeC::Ptr{nothing}
end

function __init__()
    print("\n")
    print("  |  _  __  ,__ \n")
    print("  | |_| | \\ |__ \n")
    print("\\_/ | | |_/ |__ \n")
    print("\nThis is JADE v.0.5, an interface to CASCADE,\n\n")
    print("The C-Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations!\n")
end

function acb_ode(polys::Array{acb_poly,1})
    if (length(polys) == 0)
        return
    end
    deg = 0
    order = length(polys)-1
    while iszero(polys[order+1])
        order -= 1
    end
    if order < 0
        return
    end
    A = acb_ode(order,Array{acb_poly,1}(undef,order+1),Ptr{nothing}(0))
    A.polys = polys
    finalizer(deleteC, A)
    return A
end

function setPolynomial(ode::acb_ode, index::Integer, polynomial::acb_poly)
    if ode.odeC != Ptr{nothing}(0)
        deleteC(ode)
    end
    if index > ode.order+1
        if iszero(polynomial)
            return
        end
        arr = Array{acb_poly,1}(undef,index)
        for i = 1:index-1
            if i <= ode.order+1
                arr[i] = ode.polys[i]
            else
                arr[i] = ode.polys[i].parent(0)
            end
        end
        ode.polys = arr
        ode.order = index-1
    end
    ode.polys[index] = polynomial
    # The number of polynomials might have changed
    ord = ode.order
    while iszero(ode.polys[ord+1])
        ord -= 1
    end
    if ord != ode.order
        ode.order = ord
        arr = typeof(ode.polys)(undef,ode.order+1)
        for i = 1:ode.order+1
            arr[i] = ode.polys[i]
        end
        ode.polys = arr
    end
    return
end

function translateC(ode::acb_ode)
    if ode.odeC != Ptr{nothing}(0)
        deleteC(ode)
    end
    deg = 0
    for p in ode.polys
        if deg < degree(p)
            deg = degree(p)
        end
    end
    A = ccall( (:acb_ode_setup_blank,"libcascade"), Ptr{nothing}, (Cint,Cint), deg, ode.order)
    for i = 1:ode.order+1
        ccall( (:acb_ode_set_poly,"libcascade"), Cint, (Ptr{nothing},acb_poly,Cint), A, ode.polys[i], i-1)
    end
    ode.odeC = A
    return A
end

function deleteC(ode::acb_ode)
    if ode.odeC == Ptr{nothing}(0)
        return
    end
    ccall((:acb_ode_clear,"libcascade"), Cvoid, (Ptr{nothing},), ode.odeC)
    ode.odeC = Ptr{nothing}(0)
    return
end

function setInitialValues(ode::acb_ode,poly::acb_poly)
    if ode.odeC == Ptr{nothing}(0)
        translateC(ode);
    end
    ccall((:acb_set_initial,"libcascade"), Cvoid, (Ptr{nothing},acb_poly), ode.odeC, poly)
end

function powerSeries(ode::acb_ode,target::acb)
    if ode.odeC == Ptr{nothing}(0)
        translateC(ode);
    end
    polyRing = ode.polys[1].parent
    p =  ccall((:find_power_series_julia,"libcascade"),acb_poly,(Ptr{nothing},acb,Cint),ode.odeC,target,polyRing.base_ring.prec)
    return polyRing(p)
end

function monodromy(ode::acb_ode,z0=0)
    if ode.odeC == Ptr{nothing}(0)
        translateC(ode);
    end
    Ring = ode.polys[1].parent.base_ring
    S = MatrixSpace(Ring,ode.order,ode.order)
    mono = S(1)
    ccall((:find_monodromy_matrix,"libcascade"),Cvoid,(acb_mat,Ptr{nothing},acb,Cint),mono,ode.odeC,Ring(z0),Ring.prec)
    mono.base_ring = Ring
    deleteC(ode)
    return mono
end

function acb_ode_legendre(R::AcbPolyRing,n::Integer)
    ode = acb_ode([R(n*(n+1)),R([0,-2]),R([1,0,-1])])
    ode.odeC = ccall((:acb_ode_legendre,:libcascade),Ptr{nothing},(Cint,),n)
    return ode
end

function acb_ode_bessel(R::AcbPolyRing,nu::acb)
    p = R([0,0,1])
    ode = acb_ode([p-nu*nu,R([0,1]),R([0,0,1])])
    ode.odeC = ccall((:acb_ode_bessel,:libcascade),Ptr{nothing},(acb,Cint),nu,R.base_ring.prec)
    return ode
end

function acb_ode_hypgeom(R::AcbPolyRing,a::acb,b::acb,c::acb)
    ode = acb_ode([R(-a*b),R([c,(a+b+1)]),R([0,1,-1])])
    ode.odeC = ccall((:acb_ode_hypgeom,:libcascade),Ptr{nothing},(acb,acb,acb,Cint),a,b,c,R.base_ring.prec)
    return ode
end

end
