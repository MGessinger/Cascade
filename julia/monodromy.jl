module Monodromy

using Nemo

struct acb_ode
    order::BigInt
    degree::BigInt
    series::acb_poly
    polys::Array{acb,2}
end

function acb_ode(polys::Array{acb_poly,1},initial::acb_poly,order::BigInt)
    return ccall(("acb_ode_init","monodromy"),Ptr{Nothing},(Ptr{acb_poly},acb_poly,BigInt),Ref(polys,1),initial,order)
end

function findPowerSeries(ODE::acb_ode,z::acb)
    ccall(("find_power_series","monodromy"),UInt64,(Ptr{acb_ode},acb,BigInt),Ref(ODE),z,z.parent.prec)
end

end
