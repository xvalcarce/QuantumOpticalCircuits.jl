import Base: copy
import LinearAlgebra: I,inv

mutable struct GaussianState <: AbstractState
	d::Vector{Float64}
	σ::Matrix{Float64}
end

mutable struct PseudoGaussianState <: AbstractState
	prob::Vector{Float64}
	states::Vector{GaussianState}
end

PseudoGaussianState(state::GaussianState) = PseudoGaussianState([1.0],[state])

GaussianState(modes::Int) = GaussianState(zeros(2modes),Matrix{Float64}((1/4)I,2modes,2modes))
PseudoGaussianState(modes::Int) = PseudoGaussianState([1.0],[GaussianState(modes)])

vacuum(n::Int) = GaussianState(n)
coherent(α::Float64) = GaussianState([real(α),imag(α)],Matrix{Float64}((1/4)I,2,2))

nmode(state::GaussianState) = Int(length(state.d)/2)
nmode(state::PseudoGaussianState) = Int(length(state.states[1].d)/2)

# Utilities

function copy(state::GaussianState)
	d_ = copy(state.d)
	σ_ = copy(state.σ)
	state_ = GaussianState(d_,σ_)
	return state_
end

function copy(state::PseudoGaussianState)
	prob_ = copy(state.prob)
	states_ = [copy(s) for s in state.states]
	state_ = PseudoGaussianState(prob_,states_)
	return state_
end

function M_matrix(state::GaussianState)
	m = inv(state.σ)
	return m
end

# Operations on state

Base.:(==)(s::GaussianState,t::GaussianState) = (s.d == t.d) && (s.σ == t.σ)

function Base.kron(s::GaussianState,t::GaussianState)
	ns = nmode(s)
	nt = nmode(t)
	n = ns+nt
	d = vcat(s.d,t.d)
	σ = zeros(2n,2n)
	σ[1:2ns,1:2ns] = s.σ
	base = 2ns
	σ[base+1:base+2nt,base+1:base+2nt] = t.σ
	st = GaussianState(vcat(s.d,t.d),σ)
	return st
end

function ptrace!(state::GaussianState,mode::Int)
	idx = 2mode
	deleteat!(state.d,[idx-1,idx])
	# hacky and resource heavy TODO
	σ_ = [state.σ[1:idx-2,1:idx-2] state.σ[1:idx-2,idx+1:end]
		  state.σ[idx+1:end,1:idx-2] state.σ[idx+1:end,idx+1:end]]
	state.σ = σ_
end

function ptrace!(state::GaussianState,mode::Vector{Int})
	for m in mode
		ptrace!(state,m)
		mode .-=1
	end
end
