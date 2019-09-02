import Base.+, Base.-, Base.*, Base./, Base.==

function Base.show(io::IO, A::diffEq)
    print("Order: ",A.order)
    for i = 1:length(A.polys)
        print("\npolys[$i] = ",A.polys[i])
    end
end

function ==(a::diffEq,b::diffEq)
    if (a.order != b.order)
        return false
    end
    for i = 1:a.order+1
        if (a.polys[i] != b.polys[i])
            return false
        end
    end
    return true
end

function +(a::jade_ode{P},b::jade_ode{P}) where P <: PolyElem
    x = b.order
    if (a.order > b.order)
        x = a.order
    end
    V = Vector{P}(undef,x+1)
    for i = 1:x+1
        if i > a.order+1
            V[i] = b.polys[i]
        elseif i > b.order+1
            V[i] = a.polys[i]
        else
            V[i] = a.polys[i] + b.polys[i]
        end
    end
    return jade_ode{P}(V)
end

function *(a::diffEq,n::FieldElem)
    x = a.order
    V = Vector{PolyElem}(undef,x+1)
    try
        for i = 1:x+1
            V[i] = a.polys[i]*n
        end
    catch e
        if isa(e,MethodError)
            @warn "Cannot multiply $(typeof(a)) and $(typeof(n)) (probably due to limitations in Nemo's type promotion)"
            V = Vector{fmpz_poly}([])
        else
            throw(e)
        end
    end
    return jade_ode(x,V,C_NULL)
end

function *(a::diffEq,n::Number)
    x = a.order
    V = typeof(a.polys)(undef,x+1)
    for i = eachindex(a.polys)
        V[i] = a.polys[i]*n
    end
    return jade_ode(x,V,C_NULL)
end

# Because multiplication is commutative, we introduce two methods
function *(n::Union{FieldElem,Number},a::diffEq)
    return a*n
end

function -(a::diffEq, b::diffEq)
    c = b*(-1)
    return a+c
end

function /(a::diffEq,n::Union{Number,FieldElem})
    return a*inv(n)
end
