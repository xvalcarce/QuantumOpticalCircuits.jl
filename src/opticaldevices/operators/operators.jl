import LinearAlgebra: I
import SparseArrays: SparseMatrixCSC, sparse
import QuantumOptics: create, destroy, number, identityoperator, tensor, dense, displace, dagger, SparseOperator

abstract type AbstractOperator end

abstract type SingleModeOperator <: AbstractOperator end
abstract type TwoModeOperator <: AbstractOperator end

include("./single_mode/displacement.jl")
include("./single_mode/phase_shifter.jl")
include("./single_mode/paulis.jl")
include("./single_mode/single_mode_squeezer.jl")

include("./two_modes/two_mode_squeezer.jl")
include("./two_modes/beam_splitter.jl")
include("./two_modes/swap.jl")

nmode(op::SingleModeOperator) = 1
nmode(op::TwoModeOperator) = 2

function mat1toN(mat::Matrix{Float64},N::Int,idx::NTuple{M,Int}) where {M}
	m = sparse(N,N,1.0)
	if M == 1
		idx = idx[1]
		m[2idx-1:2idx,2idx-1:2idx] = mat
	elseif M ==2 
		if idx[2]-idx[1] == 1
			[2idx[1]-1:2idx[2],2idx[1]-1:2idx[2]] = mat
		else
			#Need some hacky stuff here
		end
	else
		throw(">2 modes operators are not implemented")
	end
	return m
end

function (op::SingleModeOperator)(state::GaussianState,n::NTuple{1,Int})
	N = nmode(state)
	mat = mat(op)
	mat = mat1toN(mat,N,n)
	state.d = mat*state.d
	state.σ = mat*state.σ*mat'
end
