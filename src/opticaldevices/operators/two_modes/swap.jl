export SWAP

struct SWAP <: TwoModeOperator
end

function mat(::SWAP)
	m = sparse([4,3,2,1],
		[1,2,3,4],
		[1,-1,1,-1.])
	return m
end

function mat(::SWAP,dim::FockBasis)
	return m
end
