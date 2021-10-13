export PhaseShifter

"""
PhaseShifter

	Phase shifter operation
"""

struct PhaseShifter <: SingleModeOperator
	ϕ::Real
end

function mat(ps::PhaseShiter)
	cϕ = cos(ps.ϕ)
	sϕ = sin(ps.ϕ)
	m = [cϕ -sϕ
		 sϕ cϕ]
	return m
end

function mat(ps::PhaseShifter,dim::FockBasis)
    x = ps.ϕ*im*number(dim)
    m = exp(dense(x))
    return m
end
