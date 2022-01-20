import LinearAlgebra: I, kron
import SparseArrays: SparseMatrixCSC, sparse
import QuantumOptics: create, destroy, number, identityoperator, tensor, dense, displace, dagger, SparseOperator, Ket
import QuantumOptics: AbstractOperator as QOAbsOp

export SingleModeOperator, TwoModeOperator

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

function Ω(n::Int)
	idn = Matrix{Float64}(I,n,n)
	Ω = kron(idn,[0 1
				  -1 0])
	return Ω
end

function issymplectic(mat::Union{Matrix{Float64},SparseMatrixCSC},tol=15)
	n = Int(size(mat)[1]/2)
	M = Ω(n)
	lhs = round.(transpose(mat)*M*mat,digits=tol)
	return lhs == M
end

function mat1toN(mat::Matrix{Float64},N::Int,idx::NTuple{1,Int})
	if N == 1
		return mat
	else
		m = sparse(1.0I,2N,2N)
		idx = idx[1]
		m[2idx-1:2idx,2idx-1:2idx] = mat
		return m
	end
end

function mat1toN(mat::Union{Matrix{Float64},SparseMatrixCSC},N::Int,idx::NTuple{2,Int})
	if N == 2
		return mat
	else
		m = sparse(1.0I,2N,2N)
		if idx[2]-idx[1] == 1
			m[2idx[1]-1:2idx[2],2idx[1]-1:2idx[2]] = mat
		else
			#Need some hacky stuff here
		end
		return m
	end
end

function mat1toN(mat::QOAbsOp,N::Int,n::NTuple{1,Int})
	if N == 1
		return mat
	else
		idd = identityoperator(mat.basis_l)
		n = n[1]
		up = n == 1 ? () : (idd for i in 1:(n-1))
		down = n == N ? () : (idd for i in 1:(N-n))
		mat = tensor(up...,mat,down...)
		return mat
	end
end


function mat1toN(mat::QOAbsOp,N::Int,n::NTuple{2,Int})
	if N == 2
		return mat
	else
		idd = identityoperator(mat.basis_l.bases[1])
		Δn = abs(n[2]-n[1])
		if Δn == 1 
			n = n[1]
			up = n == 1 ? () : (idd for i in 1:(n-1))
			down = n == N ? () : (idd for i in 1:(N-n-1))
			mat = tensor(up...,mat,down...)
		else
			# TODO: find a more elegant way
			n = sort(collect(n))
			up = (idd for i in 1:(n[1]-1))
			down = (idd for i in (n[2]+1):N)
			mat = tensor(mat, (idd in 1:(n[2]-n[1]))...)
			swap_v = Vector{QOAbsOp}()
			swap_op = sparse(BS(π/2,mat.basis_l))
			for m in 0:Δn-2
				idd_swap_up = (idd for i in 1:(Δn-1-m))
				idd_swap_do = (idd for i in 1:m)
				s = tensor(idd_swap_up...,swap_op,idd_swap_do...)
				push!(swap_v,s)
			end
			mat = *(swap_v...,mat,reverse(swap_v)...)
			mat = tensor(up...,mat,down...)
		end
		return mat
	end
end

function (op::Displacement)(state::GaussianState,n::NTuple{1,Int})
	v = vec(op)
	idx = n[1]
	state.d[2idx-1] += v[1]
	state.d[2idx] += v[2]
	return state
end
	
function (op::SingleModeOperator)(state::GaussianState,n::NTuple{1,Int})
	N = nmode(state)
	m = mat(op)
	m = mat1toN(m,N,n)
	state.d = m*state.d
	state.σ = m*state.σ*m'
	return state
end

function (op::TwoModeOperator)(state::GaussianState,n::NTuple{2,Int})
	N = nmode(state)
	m = mat(op)
	m = mat1toN(m,N,n)
	state.d = m*state.d
	state.σ = m*state.σ*m'
	return state
end

function (op::AbstractOperator)(state::PseudoGaussianState,n)
	for s in state.states
		op(s,n)
	end
	return state
end

function (op::SingleModeOperator)(state::FockState,n::NTuple{1,Int})
	N = nmode(state)
	m = mat(op,state.dim)
	m = mat1toN(m,N,n)
	state.ρ = isa(state.ρ, Ket) ? m*state.ρ : m*state.ρ*dagger(m)
	return state
end

function (op::TwoModeOperator)(state::FockState,n::NTuple{2,Int})
	N = nmode(state)
	m = mat(op,state.dim)
	m = mat1toN(m,N,n)
	state.ρ = isa(state.ρ, Ket) ? m*state.ρ : m*state.ρ*dagger(m)
	return state
end
