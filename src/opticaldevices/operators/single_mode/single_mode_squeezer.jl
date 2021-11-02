export SingleModeSqueezer, SMS

""" 
Single-mode Squeezer (SMS)

"""

struct SingleModeSqueezer <: SingleModeOperator
	r::Real
	ϕ::Real
end

SMS(r::Real,ϕ::Real) = SingleModeSqueezer(r,ϕ)
SMS(r::Real) = SingleModeSqueezer(r,0.0)

function mat(sms::SingleModeSqueezer)
	shr = sinh(sms.r)
	chr = cosh(sms.r)
	sϕshr_ = -sin(sms.ϕ)*shr
	cϕshr = cos(sms.ϕ)*shr
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
