function jade_ode_legendre(R::PolyRing,n::Integer)
	ode = jade_ode([R(n*(n+1)),R([0,-2]),R([1,0,-1])])
	ode.odeC = ccall((:acb_ode_legendre,:libcascade),Ptr{Nothing},(Cint,),n)
	return ode
end

function jade_ode_bessel(R::AcbPolyRing,nu::acb)
	p = R([0,0,1])
	ode = acb_ode([p-nu*nu,R([0,1]),R([0,0,1])])
	N = Ref{acb}(nu)
	ode.odeC = ccall((:acb_ode_bessel,:libcascade),Ptr{Nothing},(Ref{acb},Cint),N,R.base_ring.prec)
	return ode
end

function jade_ode_hypgeom(R::AcbPolyRing,a::acb,b::acb,c::acb)
	ode = acb_ode([R(-a*b),R([c,(a+b+1)]),R([0,1,-1])])
	A = Ref{acb}(a)
	B = Ref{acb}(b)
	C = Ref{acb}(c)
	ode.odeC = ccall((:acb_ode_hypgeom,:libcascade),Ptr{Nothing},(Ref{acb},Ref{acb},Ref{acb},Cint),A,B,C,R.base_ring.prec)
	return ode
end

