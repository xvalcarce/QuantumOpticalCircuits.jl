export AbstractState,
	   FockState,
	   GaussianState,
	   PseudoGaussianState,
	   vacuum, coherent,
	   ptrace!,
	   nmode

import QuantumOptics: expect

abstract type AbstractState end

include("states_fock.jl")
include("states_gaussian.jl")

function GaussianState(state::FockState)
	idd = identityoperator(state.dim)
	c = create(state.dim)
	d = destroy(state.dim)
	x = (c+d)/2
	p = im*(c-d)/2
	Q = []
	for i in 1:state.n_mode
		for q in [x,p]
			vec = [idd for i in 1:state.n_mode]
			vec[i] = q
			push!(Q,tensor(vec...))
		end
	end
	ρ = isa(state.ρ,Ket) ? state.ρ ⊗ dagger(state.ρ) : state.ρ 
	d_ = [expect(ρ,q) for q in Q]
	σ_ = zeros(Complex,2state.n_mode,2state.n_mode)
	for (i,qi) in enumerate(Q)
		for (j,qj) in enumerate(Q)
			σ_[i,j] = .5*expect(ρ,qi*qj+qj*qi)-d_[i]d_[j]
		end
	end
	gs = GaussianState(d_,σ_)
	return gs
end
