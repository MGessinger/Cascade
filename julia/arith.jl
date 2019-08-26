import Base.+, Base.-, Base.*, Base./, Base.==

function Base.show(io::IO, A::acb_ode)
    if A.odeC != C_NULL
        ccall((:acb_ode_dump,:libcascade),Cvoid,(Ptr{Nothing},Ptr{Nothing}),A.odeC,C_NULL)
        return
    end
    print("Order: ",A.order)
    for i = 1:length(A.polys)
        print("\npolys[",i,"] = ",A.polys[i])
    end
end

function ==(a::acb_ode,b::acb_ode)
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

function +(a::acb_ode,b::acb_ode)
    x = a.order
    if (b.order > x)
       x = b.order
    end
    V = Vector{acb_poly}(undef,x+1)
    for i = 1:x+1
       if (i <= a.order+1 && i <= b.order+1)
           V[i] = a.polys[i] + b.polys[i]
       elseif (i > b.order+1)
           V[i] = a.polys[i]
       else
           V[i] = b.polys[i]
       end
    end
    return acb_ode(V)
end

function *(a::acb_ode,n::Union{FieldElem,Number})
    x = a.order
    V = Vector{acb_poly}(undef,x+1)
    for i = 1:x+1
        V[i] = a.polys[i]*n
    end
    return acb_ode(V)
end

# Because multiplication is commutative, we introduce two methods
function *(n::Union{FieldElem,Number},a::acb_ode)
	return a*n
end

function -(a::acb_ode, b::acb_ode)
    c = b*(-1)
    return a+c
end

function /(a::acb_ode,n::Union{Number,FieldElem})
	return a*(1/n)
end
