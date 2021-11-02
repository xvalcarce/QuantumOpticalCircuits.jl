export SWAP, swap

struct SWAP <: TwoModeOperator
end

swap() = SWAP()

function mat(::SWAP)
	m = sparse([4,3,2,1],
		[1,2,3,4],
		[1,-1,1,-1.])
	return m
end

function mat(::SWAP,dim::FockBasis)	
	x = (im*π/2)*tensor(create(dim),destroy(dim))+tensor(destroy(dim),create(dim))
	m = exp(dense(x))
	return m
end
