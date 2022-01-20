export BeamSplitter, BS

"""
Beam splitter

"""

struct BeamSplitter <: TwoModeOperator
	θ::Real
end

BS(θ::Real) = BeamSplitter(θ)

"""
	BeamSplitter(θ)

	Return a beam splitter operation with transmitivity cos(θ)
"""

function mat(bs::BeamSplitter)
	cθ = cos(bs.θ)
	sθ = sin(bs.θ)
	r = sparse([1,2,3,4,4,3,2,1],
		[1,2,3,4,1,2,3,4],
		[cθ,cθ,cθ,cθ,sθ,-sθ,sθ,-sθ])
	return r
end

function mat(bs::BeamSplitter,dim::FockBasis)
    x = tensor(create(dim),destroy(dim))+tensor(destroy(dim),create(dim))
    x *= im*bs.θ
    m = exp(dense(x))
    return m
end

Base.:(==)(lhs::BeamSplitter,rhs::BeamSplitter) = lhs.θ == rhs.θ
