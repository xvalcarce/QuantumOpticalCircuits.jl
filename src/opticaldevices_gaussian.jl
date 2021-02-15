import SparseArrays: SparseMatrixCSC, sparse

# Utils

id2 = [1 0
	  0 1]

function s_matrix(θ::Float64)
	r = [cos(θ) sin(θ)
		sin(θ) -cos(θ)]
	return r
end

# Single-mode optical devices

# Phase-shifter
function ps(ϕ::Float64)
	r = [cos(ϕ) -sin(ϕ)
		 sin(ϕ) cos(ϕ)]
	return r
end

# Single-mode squeezer
function sms_re(g::Float64)
	r = [exp(-g) 0
		 0 exp(g)]
	return r
end

function sms_im(g::Float64)
	r = [cosh(g) -sinh(g)
		 -sinh(g) cosh(g)]
	return r
end

function sms(g::Complex)
	rg = real(g)
	ig = imag(g)
	r = sms_re(rg)*(sms_im(ig))
	return r
end

function sms(r::Float64,ϕ::Float64)
	sϕ = s_matrix(ϕ)
	r = cosh(r)*id2-sinh(r)*sϕ
	return r
end

# Displacement
function disp_re(α::Float64,d::Vector{Complex{Float64}},mode::Int)
	d_ = d
	d_[2mode-1] += α
	return d_
end

function disp_im(α::Float64,d::Vector{Complex{Float64}},mode::Int)
	d_ = d
	d_[2mode] += α
	return d_
end

function disp(α::Complex,d::Vector{Complex{Float64}},mode::Int)
	d_ = d
	d_[2mode-1] += real(α)
	d_[2mode] += imag(α)
	return d_
end

# Two-mode optical devices

# Two-mode squeezer
function tms(r::Float64,ϕ::Float64)
	chr = cosh(r)
	shr = sinh(r)
	sϕ = s_matrix(ϕ)
	shrϕ = shr*sϕ
	chrϕ = chr*id2
	r = [chrϕ shrϕ
		shrϕ chrϕ]	
	return r
end

# Beam-splitter
# TODO: Implement phase angle
function bs(θ::Float64)
	cθ = cos(θ)
	sθ = sin(θ)
	r = [cθ 0 0 -sθ
		0 cθ sθ 0
		0 -sθ cθ 0
		sθ 0 0 cθ]
	return r
end

# Swap gate
function swap()
	r = sparse([4,3,2,1],
			[1,2,3,4],
			[1,-1,1,-1.])
	return r
end

# Application of gate on Gaussian state

function apply(od::OpticalDevice,state::GaussianState)
	state_ = copy(state)
	apply!(od,state_)
	return state_
end

function apply!(od::OpticalDevice,state::GaussianState)
	if od.optdev ∈ [disp,disp_re,disp_im]
		state.d = od.optdev(od.param,state.d,od.mode)
	else	
		dg = collect(1:length(state.d))
		r = sparse(dg,dg,1.0)
		r_ = od.optdev(od.param...)
		if typeof(od.mode) == Int
			r[2od.mode-1:2od.mode,2od.mode-1:2od.mode] = r_
		elseif length(od.mode) == 2
			Δmode = abs(od.mode[2]-od.mode[1])
			if Δmode == 1
				r[2od.mode[1]-1:2od.mode[2],2od.mode[1]-1:2od.mode[2]] = r_
			else 
				# Ensuring full connectivity
				# TODO: possible speed up swap mode (n,m) directly?
				swaps = Vector{SparseMatrixCSC{Float64,Int64}}()
				swap_ = swap()
				dom,upm = od.mode[2]>od.mode[1] ? od.mode : sort(od.mode) 
				for m in 1:Δmode-1
					m_ = upm-m
					s = sparse(dg,dg,1.0)
					s[2m_-1:2(m_+1),2m_-1:2(m_+1)] = swap_
					push!(swaps,s)
				end
				r[2dom-1:2(dom+1),2dom-1:2(dom+1)] = r_
				r = *(swaps...,r,reverse(swaps)...)
			end
		else
			throw(">2 modes operations are not implemented")
		end
		state.d = r*state.d
		state.σ = r*state.σ*r'
	end
end

function apply(od::OpticalDevice,state::PseudoGaussianState)
	state_ = copy(state)
	for (idx,s) in state_.states
		apply!(od,s)
	end
	return state_
end

function apply!(od::OpticalDevice,state::PseudoGaussianState)
	for s in state.states
		apply!(od,s)
	end
end
