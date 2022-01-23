export OpticalDevice

include("./operators/operators.jl")

struct OpticalDevice{M,OP<:AbstractOperator}
	operator::OP
	mode::NTuple{M,Int}

	function OpticalDevice(operator::OP,mode::NTuple{M,Int}) where {M,OP<:AbstractOperator}
		@assert nmode(operator) == M "Operator acts on $(nmode(operator)) mode(s), $M were given."
		return new{M,OP}(operator,mode)
	end
end

(op::AbstractOperator)(i::Int) = OpticalDevice(op,(i,))
(op::AbstractOperator)(i::Int,j::Int) = OpticalDevice(op,(i,j))
(op::AbstractOperator)(m::Vector{Int}) = OpticalDevice(op,Tuple(m))
(op::AbstractOperator)(m::NTuple{M,Int} where {M}) = OpticalDevice(op,m)

function (od::OpticalDevice)(state::AbstractState)
	@assert all(od.mode .≤ nmode(state)) "OpticalDevice acts on a mode larger than the state"
	return od.operator(state,od.mode)
end
