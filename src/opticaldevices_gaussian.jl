function ps(ϕ::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	r = [cos(ϕ) sin(ϕ)
		 -sin(ϕ) cos(ϕ)]
	d_ = r*d
	σ_ = r*σ*r'
	return d_,σ_
end

function sms_re(g::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	r = [exp(-g) 0
		 0 exp(g)]
	d_ = r*d
	σ_ = r*σ*r'
	return d_,σ_
end

function sms_im(g::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	r = [cosh(g) -sinh(g)
		 -sinh(g) cosh(g)]
	d_ = r*d
	σ_ = r*σ*r'
	return d_,σ_
end

function sms(g::Complex,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	d_,σ_ = sms_re(real(g),d,σ)
	d_,σ_ = sms_im(imag(g),d_,σ_)
end

function disp_re(α::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	d_ = d
	d_[1] += α
	σ_ = σ
	return d_,σ_
end

function disp_im(α::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	d_ = d
	d_[2] += α
	σ_ = σ
	return d_,σ_
end

function disp(α::Complex,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	d_,σ_ = disp_re(real(α),d,σ)
	d_,σ_ = disp_im(imag(α),d_,σ)
	return d_,σ_
end

function tms(θ::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	r = [cosh(θ) 0 sinh(θ) 0
		 0 cosh(θ) 0 sinh(θ)
		 sinh(θ) 0 cosh(θ) 0
		 0 sinh(θ) 0 cosh(θ)]
	d_ = r*d
	σ_ = r*σ*r'
	return d_,σ_
end

function bs(θ::Float64,d::Array{Complex{Float64},1},σ::Array{Complex{Float64},2})
	r = [cos(θ) 0 sin(θ)*im 0
		 0 cos(θ) 0 sin(θ)*im
		 sin(θ)*im 0 cos(θ) 0
		 0 sin(θ)*im 0 cos(θ)]
	d_ = r*d
	σ_ = r*σ*r'
	return d_,σ_
end

function apply(od::OpticalDevice,state::GaussianState)
	state_ = copy(state)
	if typeof(od.mode) == Int
		d,σ = extract_one_mode(state_,od.mode)
		d_,σ_ = od.optdev(od.param,d,σ)
		update!(state_,od.mode,d_,σ_)
	else
		d,σ = extract_two_modes(state_,od.mode)
		d_,σ_ = od.optdev(od.param,d,σ)
		update!(state_,od.mode,d_,σ_)
	end
	return state_
end

function apply!(od::OpticalDevice,state::GaussianState)
	if typeof(od.mode) == Int
		d,σ = extract_one_mode(state,od.mode)
		d_,σ_ = od.optdev(od.param,d,σ)
		update!(state,od.mode,d_,σ_)
	else
		d,σ = extract_two_modes(state,od.mode)
		d_,σ_ = od.optdev(od.param,d,σ)
		update!(state,od.mode,d_,σ_)
	end
end
