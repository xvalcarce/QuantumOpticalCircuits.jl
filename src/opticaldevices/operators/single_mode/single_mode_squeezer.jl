export SingleModeSqueezer

""" 
Single-mode Squeezer (SMS)

"""

struct SingleModeSqueezer <: SingleModeOperator
	r::Real
	ϕ::Real
end

function mat(sms::SingleModeSqueezer)
	shr = sinh(r)
	chr = cosh(r)
	sϕshr_ = -sin(ϕ)*shr
	cϕshr = cos(ϕ)*shr
	m = [chr-cϕshr sϕshr_
		 sϕshr_ chr+cϕshr]
	return m
end

function mat(sms::SingleModeSqueezer,dim::FockBasis)
	g = sms.r*exp(sms.ϕ*im)
	x = .5*(conj(g)*(destroy(dim)^2)-g*(create(dim)^2))
    m = exp(dense(x))
    return m
end
