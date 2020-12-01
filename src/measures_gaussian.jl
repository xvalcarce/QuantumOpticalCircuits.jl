import LinearAlgebra: I,det,inv

function p_noclick(state::GaussianState,mode::Int,η::Float64)
	d = state.d
	M = M_matrix(state)
	F = zeros(size(M)...)
	F[2mode-1:2mode,2mode-1:2mode] = (4*(1-η)/(1+η))*Matrix{Float64}(I,2,2)
	p_nc = (2*√(det(M)))/((1+η)*√(det(M+F)))
	p_nc *= exp(-.5*d'*(M-M*inv(M+F)*M)*d)
	p_nc = real(p_nc)
	return p_nc
end

function p_noclick!(state::GaussianState,mode::Int,η::Float64)
	d = state.d
	M = M_matrix(state)
	F = zeros(size(M)...)
	F[2mode-1:2mode,2mode-1:2mode] = (4*(1-η)/(1+η))*Matrix{Float64}(I,2,2)
	p_nc = (2*√(det(M)))/((1+η)*√(det(M+F)))
	p_nc *= exp(-.5*d'*M*d+.5*d'*M*inv(M+F)*M*d)
	p_nc = real(p_nc)
	σ_inv = inv(state.σ)
	d_ = inv(σ_inv+F)*σ_inv*state.d
	σ_ = inv(σ_inv+F)
	state.d = d_
	state.σ = σ_
	trace_mode!(state,mode)
	return p_nc
end

function correlator(state::GaussianState,i::Int,j::Int,η::Float64)
	if i==j
		throw(ArgumentError("Modes i,j need to be different"))
	elseif i>j
		i,j = j,i
	end
	state_ = copy(state)
	p_b = p_noclick(state,j,η)
    p_a = p_noclick!(state_,i,η)
    p_b_cond = p_noclick(state_,j-1,η)
	corr = 1 - 2*(p_a+p_b) + 4*p_a*p_b_cond
    return corr
end
