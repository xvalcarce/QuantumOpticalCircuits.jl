export TwoModeSqueezer

struct TwoModeSqueezer <: TwoModeOperator
	r::Real
	ϕ::Real
end

function mat(tms::TwoModeSqueezer)
	cϕ = cos(tms.ϕ)
	sϕ = sin(tms.ϕ)
	shr = sinh(tms.r)
	diag = cosh(tms.r)*[1 0
					  0 1]
	phase = shr*[cϕ sϕ
				 sϕ -cϕ]
	m = [diag phase
		 phase diag]
	return m
end

function mat(tms::TwoModeSqueezer,dim::FockBasis)
	z = tms.r*exp(tms.ϕ*im)
	x = z*tensor(create(dim),create(dim))-conj(z)*tensor(destroy(dim),destroy(dim))
    m = exp(dense(x))
	return m
end
