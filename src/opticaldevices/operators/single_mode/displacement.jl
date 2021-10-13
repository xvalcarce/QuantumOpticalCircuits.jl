export Displacement

struct Displacement <: SingleModeOperator
	α::Complex
end	

Displacement(α::Real) = Displacement(Complex(α))

function mat(d::Displacement,dim::FockBasis)
	m = displace(dim,d.α)
	return m
end

function mat(d::Displacement)
	@warning "Displacement in gaussian representation can not be seen as a matrix operation.\nReturning identity as a fallback"
	return Matrix{1.0I, 2, 2}
end

function vec(d::Displacement)
	v = [real(d.α)
		 imag(d.α)]
	return v
end
