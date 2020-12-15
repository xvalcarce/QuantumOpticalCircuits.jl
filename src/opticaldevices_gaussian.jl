import SparseArrays: SparseMatrixCSC, sparse

function ps(ϕ::Float64)
	r = [cos(ϕ) sin(ϕ)
		 -sin(ϕ) cos(ϕ)]
	return r
end

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
	d_ = disp_re(real(α),d,mode)
	d_ = disp_im(imag(α),d_,mode)
	return d_
end

function tms(θ::Float64)
	shθ = sinh(θ)
	chθ = cosh(θ)
	r = [chθ 0 shθ 0
		 0 chθ 0 -shθ
		 shθ 0 chθ 0
		 0 -shθ 0 chθ]
	return r
end

function bs(θ::Float64)
	cθ = cos(θ)
	sθ = sin(θ)
	r = [cθ 0 0 -sθ
		0 cθ sθ 0
		0 -sθ cθ 0
		sθ 0 0 cθ]
	return r
end

function swap()
	r = sparse([4,3,2,1],
			[1,2,3,4],
			[1,-1,1,-1.])
	return r
end

function apply(od::OpticalDevice,state::GaussianState)
	if od.optdev ∈ [disp,disp_re,disp_im]
		d = od.optdev(od.param,state.d,od.mode)	
		σ = state.σ
	else
		dg = collect(1:length(state.d))
		r = sparse(dg,dg,1.0)
		r_ = od.optdev(od.param)
		if typeof(od.mode) == Int
			r[2od.mode-1:2od.mode,2od.mode-1:2od.mode] = r_
		elseif length(od.mode) == 2
			if od.mode[2]-od.mode[1]==1
				r[2od.mode[1]-1:2od.mode[2],2od.mode[1]-1:2od.mode[2]] = r_
			else
				# Ensuring full connectivity
				for (idx,i) in enumerate(od.mode)
					for (jdx,j) in enumerate(od.mode)
						r[2i-1:2i,2j-1:2j] = r_[2idx-1:2idx,2jdx-1:2jdx]
					end
				end
			end
		else
			throw("2+ mode operation are not implemented")
		end
		d = r*state.d
		σ = r*state.σ*r'
	end
	state_ = GaussianState(d,σ)
	return state_
end

function apply!(od::OpticalDevice,state::GaussianState)
	if od.optdev ∈ [disp,disp_re,disp_im]
		state.d = od.optdev(od.param,state.d,od.mode)
	else	
		dg = collect(1:length(state.d))
		r = sparse(dg,dg,1.0)
		r_ = od.optdev(od.param)
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
