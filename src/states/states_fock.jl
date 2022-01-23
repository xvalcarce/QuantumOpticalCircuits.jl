import Base: copy
import QuantumOptics: FockBasis, coherentstate, tensor, Ket, Operator, ptrace

mutable struct FockState <: AbstractState
	n_mode::Int
	cutoff_dim::Int
	dim::FockBasis
	ρ::Union{Operator,Ket}
end

function FockState(n_mode::Int,cutoff_dim::Int)
	dim = FockBasis(cutoff_dim)
	vacuum_mode = coherentstate(dim,0)
	vacuum = tensor((vacuum_mode for i in 1:n_mode)...)
	return FockState(n_mode,cutoff_dim,dim,vacuum)
end

vacuum(n::Int,d::Int) = FockState(n,d)

function coherent(α::Number,d::Int)
	dim = FockBasis(d)
	return FockState(1,d,dim,coherentstate(dim,α))
end

nmode(state::FockState) = state.n_mode

# Utilities

function copy(state::FockState)
	ρ_ = copy(state.ρ)
	state_ = FockState(state.n_mode,state.cutoff_dim)
	state_.ρ = ρ_
	return state_
end

Base.:(==)(s::FockState,t::FockState) = s.ρ == t.ρ

function Base.kron(s::FockState,t::FockState)
	n_mode = nmode(s)+nmode(t)
	@assert s.cutoff_dim == t.cutoff_dim "States need to have the same Hilbert space cutoff"
	ρ = s.ρ ⊗ t.ρ
	return FockState(n_mode,s.cutoff_dim,s.dim,ρ)
end

function ptrace!(state::FockState,mode::Int)
	state.ρ = ptrace(state.ρ,mode)
	state.n_mode -= 1
end

function ptrace!(state::FockState,mode::Vector{Int})
	state.ρ = ptrace(state.ρ,mode)
	state.n_mode -= length(mode)
end
