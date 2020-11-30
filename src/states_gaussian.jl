import Base: copy
import LinearAlgebra: I,inv

mutable struct GaussianState <: State
	d::Array{Complex{Float64},1}
	σ::Array{Complex{Float64},2}
end

vacuum_gaussian(modes::Int) = GaussianState(complex(zeros(2modes)),Matrix{Complex{Float64}}(I,2modes,2modes))

# Utilities

function copy(state::GaussianState)
	d_ = copy(state.d)
	σ_ = copy(state.σ)
	state_ = GaussianState(d_,σ_)
	return state_
end

function M_matrix(state::GaussianState)
	m = 4inv(state.σ)
	return m
end

# Operations on state

function trace_mode!(state::GaussianState,mode::Int)
	idx = 2mode
	deleteat!(state.d,[idx-1,idx])
	# hacky and resource heavy TODO
	σ_ = [state.σ[1:idx-2,1:idx-2] state.σ[1:idx-2,idx+1:end]
		  state.σ[idx+1:end,1:idx-2] state.σ[idx+1:end,idx+1:end]]
	state.σ = σ_
end

# Extract functions 1 and 2 mode decoupled from n -> for performance, we avoid to create enumerator this way

function extract_one_mode(state::GaussianState,mode::Int)
	idx = 2*mode
	d = state.d[idx-1:idx]
	σ = state.σ[idx-1:idx,idx-1:idx]
	return d,σ
end

function extract_two_modes(state::GaussianState,mode::Array{Int,1})
	idx1 = 2mode[1]
	idx2 = 2mode[2]
	d = [state.d[idx1-1:idx1]
		 state.d[idx2-1:idx2]]
	ii = state.σ[idx1-1:idx1,idx1-1:idx1]
	ij = state.σ[idx1-1:idx1,idx2-1:idx2]
	ji = state.σ[idx2-1:idx2,idx1-1:idx1]
	jj = state.σ[idx2-1:idx2,idx2-1:idx2]
	σ = [ii ij
		 ji jj]
	return d,σ
end

function extract_n_modes(state::GaussianState,mode::Array{Int,1})
	idxs = [2m for m in mode]
	l = 2length(mode)
	d = zeros(l)
	σ = zeros(l,l)
	for (i,idx1) in enumerate(idxs)
		d[2i-1:2i] = state.d[idx1-1:idx1]
		for (j,idx2) in enumerate(idxs)
			σ[2i-1:2i,2j-1:2j] = state.σ[idx1-1:idx1,idx2-1:idx2]
		end
	end
	return d,σ
end

function extract(state::GaussianState,mode::Union{Int,Array{Int,1}})
	if typeof(mode) == Int
		d,σ = extract_one_mode(state,mode)
	elseif length(mode) == 2
		d,σ = extract_two_modes(state,mode)
	else
		d,σ = extract_n_modes(state,mode)
	end
	return d,σ
end

# Update functions 1 and 2 mode decoupled from n -> for performance, we avoid to create enumerator this way

function update_one_mode!(state::GaussianState,mode::Int,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	idx = 2*mode
	state.d[idx-1:idx] = d
	state.σ[idx-1:idx,idx-1:idx] = σ
end

function update_two_modes!(state::GaussianState,mode::Array{Int,1},d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	idx1 = 2mode[1]
	idx2 = 2mode[2]
	state.d[idx1-1:idx1] = d[1:2]
	state.d[idx2-1:idx2] = d[3:4]
	state.σ[idx1-1:idx1,idx1-1:idx1] = σ[1:2,1:2]
	state.σ[idx1-1:idx1,idx2-1:idx2] = σ[1:2,3:4]
	state.σ[idx2-1:idx2,idx1-1:idx1] = σ[3:4,1:2]
	state.σ[idx2-1:idx2,idx2-1:idx2] = σ[3:4,3:4]
end

function update_n_modes!(state::GaussianState,mode::Array{Int,1},d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	idxs = [2m for m in mode]
	for (i,idx1) in enumerate(idxs)
		state.d[idx1-1:idx1] = d[2*i-1:2*i]
		for (j,idx2) in enumerate(idxs)
			state.σ[idx1-1:idx1,idx2-1:idx2] = σ[2i-1:2i,2j-1:2j]
		end
	end
end

function update!(state::GaussianState,mode::Union{Int,Array{Int,1}},d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	if typeof(mode) == Int
		update_one_mode!(state,mode,d,σ)
	elseif length(mode) == 2
		update_two_modes!(state,mode,d,σ)
	else
		update_n_modes!(state,mode,d,σ)
	end
end
