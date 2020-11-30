import QuantumOptics: FockBasis, coherentstate, tensor, Ket, Operator

mutable struct FockState <: State
	n_mode::Int
	cutoff_dim::Int
	dim::FockBasis
	ρ::Union{Operator,Ket}
	function FockState(n_mode,cutoff_dim)
		dim = FockBasis(cutoff_dim)
		vacuum_mode = coherentstate(dim,0)
		vacuum = tensor((vacuum_mode for i in 1:n_mode)...)
		new(n_mode,cutoff_dim,dim,vacuum)
	end
end

vacuum_fock(n_mode::Int,cutoff_dim::Int) = FockState(n_mode,cutoff_dim)

function copy(state::FockState)
	ρ_ = copy(state.ρ)
	state_ = FockState(state.n_mode,state.cutoff_dim)
	state_.ρ = ρ_
	return state_
end
