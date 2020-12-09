import Base: copy
import QuantumOptics: FockBasis, coherentstate, tensor, Ket, Operator, ptrace

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

function ptrace!(state::FockState,mode::Int)
	state.ρ = ptrace(state.ρ,mode)
	state.n_mode -= 1
end

function ptrace!(state::FockState,mode::Vector{Int})
	state.ρ = ptrace(state.ρ,mode)
	state.n_mode -= length(mode)
end
