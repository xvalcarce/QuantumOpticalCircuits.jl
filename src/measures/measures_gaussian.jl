import LinearAlgebra: Symmetric,I,det,inv,normalize,eigvals
import SparseArrays: spzeros
import PDMats: PDMat,det,inv

# Gaussian state

function p_noclick(state::GaussianState,mode::Int,η::Float64)
	d = state.d
	#σ = PDMat(Symmetric(state.σ))
    σ = state.σ
	M = inv(σ)
	F = spzeros(size(M)...)
	F[2mode-1:2mode,2mode-1:2mode] = (4η/(2-η))*Matrix{Float64}(I,2,2)
	M_F = PDMat(M+F)
	p_nc = (2*√(det(M)))/((2-η)*√(det(M_F)))
	p_nc *= exp(-.5*d'*(M-M*inv(M_F)*M)*d)
	p_nc = real(p_nc)
	return p_nc
end

function p_noclick!(state::GaussianState,mode::Int,η::Float64)
	d = state.d
	σ = PDMat(Symmetric(state.σ))
	M = inv(σ)
	F = spzeros(size(M)...)
	F[2mode-1:2mode,2mode-1:2mode] = (4η/(2-η))*Matrix{Float64}(I,2,2)
    M_F = PDMat(M+F)
	imf = inv(M_F)
	p_nc = 2/(2-η)
	p_nc /= √(det(σ)*det(M_F))
	p_nc *= exp(-.5*d'*(M-M*imf*M)*d)
	p_nc = real(p_nc)
	d_ = imf*M*d
	σ_ = imf
	state.d = d_
	state.σ = σ_
	ptrace!(state,mode)
	return p_nc
end

function herald_click!(state::GaussianState,mode::Int,η::Float64)
	throw(ArgumentError("Can not herald a GaussianState: the post-selected state is not Gaussian.\n
						Try with a PseudoGaussianState or a PseudoGaussianCircuit"))
end

herald_noclick!(state::GaussianState,mode::Int,η::Float64) = herald_click!(state,mode,η)

function o_(l::Vector{Int},η::Float64)
	η_factor = (4η)/(2-η)
	n = length(l)
	o = Matrix{Int}(I,2n,2n)
	for i in 1:n
		if l[i] == 1
			continue
		else
			o[2i-1,2i-1] = 0
			o[2i,2i] = 0
		end
	end
	o *= η_factor
	return o
end

function O_(state::GaussianState,l::Vector{Int},η::Float64)
	η_factor = (2/(2-η))^sum(l)
	o = o_(l,η)
	σ = PDMat(Symmetric(state.σ))
	M = inv(σ)
	d = state.d
	O = exp(-.5*d'*(M-M*inv(M+o)*M)*d)
	O /= √(det(σ)*det(M+o))
	O *= η_factor
	return real(O)
end

function E(state::GaussianState,k::Vector{Int},η::Float64)
	"""
	k -> 0 no click , 1 click
	"""
	k_0s = findall(x->x==0,k)
	lk = length(k)
	n_bitstring = lk-length(k_0s)
	if n_bitstring == 0
		S_k = [ones(Int,lk)]
	else
		S_k = collect.(Iterators.product(fill(0:1,n_bitstring)...))[:]
		for l in S_k
			for i in k_0s
				insert!(l,i,1)
			end
		end
	end
	E = 0
	for l in S_k
		E += ((-1)^(k'l))*O_(state,l,η)
	end
	return E
end

function correlator(state::GaussianState,i::Int,j::Int,η::Float64)
	if i==j
		throw(ArgumentError("Modes i,j need to be different"))
	elseif i>j
		i,j = j,i
	end
	state_ = copy(state)
	n_mode = length(state_.d)/2
	if n_mode > 2
		m = filter(e->e∉[i,j],collect(1:n_mode))
		ptrace!(state_,m)
	end
	s = copy(state_)
	p_b = p_noclick(state_,2,η)
    p_a = p_noclick!(s,1,η)
    p_b_cond = p_noclick(s,1,η)
	corr = 1 - 2*(p_a+p_b) + 4*p_a*p_b_cond
    return corr
end

function wigner(state::GaussianState,α::Vector{Float64})
	σ = PDMat(Symmetric(state.σ))
	Δα = α-state.d
	W = exp(-.5*invquad(σ,Δα))/√(det(2π*σ))
	return W
end

function wigner(state::GaussianState,αs::Vector{Vector{Float64}})
	σ = PDMat(Symmetric(state.σ))
	den = √(det(2π*σ))
	N = length(αs)
	W = zeros(N)
	Threads.@threads for i in 1:N
		Δα = αs[i]-state.d
		W[i] = exp(-.5*invquad(σ,Δα))/den
	end
	return W
end

# Pseudo-gaussian state measures

function E(state::PseudoGaussianState,k::Vector{Int},η::Float64)
	E_ = 0
	for (p,s) in zip(state.prob,state.states)
		E_ += p*E(s,k,η)
	end
	E_ *= state.norm
	return E_
end

function correlator(state::PseudoGaussianState,i::Int,j::Int,η::Float64)
	Es = [[E(state,[a,b],η) for b in 0:1] for a in 0:1]
	corr = sum([(-1)^(Int(a!=b))*Es[a][b] for a in 1:2 for b in 1:2])
    return corr
end

function p_noclick(state::PseudoGaussianState,mode::Int,η::Float64)
	p_nc = 0.0
	for (c_i,s) in zip(state.prob,state.states)
		p_nc += c_i*p_noclick(s,mode,η)
	end
	p_nc *= state.norm
	return p_nc
end

function p_noclick!(state::PseudoGaussianState,mode::Int,η::Float64)
	p_nc = 0.0
	for (idx,s) in enumerate(state.states)
		p_nc += state.prob[idx]*p_noclick!(s,mode,η)
	end
	p_nc *= state.norm
	return p_nc
end

function herald_click!(state::PseudoGaussianState,mode::Int,η::Float64,tol::Float64)
	p_c = 0
	states = Vector{GaussianState}()
	prob = Vector{Float64}()
	n_composite_state = length(state.prob)
	p_noclick(state,mode,η)
	state.norm *= 1/(1-p_noclick(state,mode,η))
	for i in 1:n_composite_state
		state_□ = copy(state.states[i])
		ptrace!(state_□,mode)
		p_nc = p_noclick!(state.states[i],mode,η)
		p_c_i = 1-p_nc
        @assert p_c_i > tol "Click probability below tolerance ($tol): $p_c_i"
		p_c += p_c_i
		append!(states,[state_□,state.states[i]])
		append!(prob,state.prob[i]*[1.0,-p_nc])
	end
	state.prob = prob
	state.states = states
	return p_c
end

herald_noclick!(state::PseudoGaussianState,mode::Int,η::Float64) = p_noclick!(state,mode,η)
