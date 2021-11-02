import Base: copy
import LinearAlgebra: I,inv

mutable struct GaussianState <: AbstractState
	d::Array{Complex{Float64},1}
	σ::Array{Complex{Float64},2}
end

mutable struct PseudoGaussianState <: AbstractState
	prob::Vector{Float64}
	states::Vector{GaussianState}
end

PseudoGaussianState(state::GaussianState) = PseudoGaussianState([1.0],[state])

GaussianState(modes::Int) = GaussianState(complex(zeros(2modes)),(1/4)Matrix{Complex{Float64}}(I,2modes,2modes))
PseudoGaussianState(modes::Int) = PseudoGaussianState([1.0],[GaussianState(modes)])

nmode(state::GaussianState) = Int(length(state.d)/2)

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
