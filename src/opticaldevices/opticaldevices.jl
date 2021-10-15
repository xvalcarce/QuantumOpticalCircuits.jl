struct OpticalDevice{M,OP<:AbstractOperator}
	operator::OP
	mode::NTuple{M,Int}

	function OpticalDevice(operator::OP,mode::NTuple{M,Int}) where {M,OP<:AbstractOperator}
		@assert nmode(OP) == M "Operator acts on $(nmode(OP)), $M were given."
		return new{M,OP}(operator,mode)
	end
end

(op::AbstractOperator)(i::Int) = OpticalDevice(op,(i,))
(op::AbstractOperator)(i::Int,j::Int) = OpticalDevice(op,(i,j))
(op::AbstractOperator)(m::Vector{Int}) = OpticalDevice(op,Tuple(m))
(op::AbstractOperator)(m::NTuple{M,Int} where {M}) = OpticalDevice(op,m)

function (od::OpticalDevice)(state::AbstractState)
	@assert any(od.mode .≤ nmode(state)) "OpticalDevice acts on a mode larger than the state"
	return od.op(state,od.mode)
end
