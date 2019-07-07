module Jade

using Nemo

global const Pol = PolyElem{acb}

export acb_ode,monodromy,powerSeries,translateC,deleteC,setPolynomial,setInitialValues

mutable struct acb_ode 
    degree::Integer
    order::Integer
    polys::Array{N,1} where N <: Pol
    odeC::Ptr{nothing}
end

function __init__()
    print("\n")
    print("  |  _  __  ,__ \n")
    print("  | |_| | \\ |__ \n")
    print("\\_/ | | |_/ |__ \n")
    print("\nThis is JADE v.0.3, an interface to CASCADE,\n\n")
    print("The C-Library for Approximative Solutions to Complex Arbitrary Precision Differential Equations!\n")
end

function acb_ode(polys::Array{T,1} where T <: Pol)
    if (length(polys) == 0)
        return
    end
    deg = 0
    order = length(polys)-1
    while iszero(polys[order+1])
        order -= 1
    end
    for p in polys
        if deg < degree(p)
            deg = degree(p)
        end
    end
    if order < 0 || deg < 0
        return
    end
    A = acb_ode(deg,order,Array{typeof(polys[1]),1}(undef,order+1),0)
    for i = 1:A.order+1
        A.polys[i] = polys[i]
    end
    return A
end

function setPolynomial(ode::acb_ode, index::Integer, polynomial::N where N <: Pol)
    if ode.odeC != Ptr{nothing}(0)
        deleteC(ode)
    end
    if index > ode.order+1
        if iszero(polynomial)
            return
        end
        arr = typeof(ode.polys)(undef,index)
        for i = 1:index-1
            if i <= ode.order+1
                arr[i] = ode.polys[i]
            else
                arr[i] = ode.polys[1].parent(0)
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
    # The value for ode.degree might have changed
    deg = 0
    for p in ode.polys
        if deg < degree(p)
            deg = degree(p)
        end
    end
    ode.degree = deg 
    return
end

function translateC(ode::acb_ode)
    if ode.odeC != Ptr{nothing}(0)
        deleteC(ode)
    end
    print("[INFO] Translating ODE to C-style struct.\n")
    A = ccall( (:acb_ode_setup_blank,"libcascade"), Ptr{nothing}, (Cint,Cint), ode.degree, ode.order)
    for i = 1:ode.order+1
        ccall( (:acb_ode_set_poly,"libcascade"), Cint, (Ptr{nothing},acb_poly,Cint), A, ode.polys[i], i-1)
    end
    ode.odeC = A
    return A
end

function deleteC(ode::acb_ode)
    ccall((:acb_ode_clear,"libcascade"), Cvoid, (Ptr{nothing},), ode.odeC)
    ode.odeC = Ptr{nothing}(0)
    return
end

function setInitialValues(ode::acb_ode,poly::N where N <: Pol)
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

function monodromy(ode::acb_ode)
    if ode.odeC == Ptr{nothing}(0)
        translateC(ode);
    end
    Ring = ode.polys[1].parent.base_ring
    S = MatrixSpace(Ring,ode.order,ode.order)
    mono = S(0)
    ccall((:find_monodromy_matrix,"libcascade"),Cvoid,(Ptr{nothing},acb_mat,Cint),ode.odeC,mono,Ring.prec)
    mono.base_ring = Ring
    deleteC(ode)
    return mono
end

end
