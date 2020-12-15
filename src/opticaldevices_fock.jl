import QuantumOptics: create, destroy, identityoperator, tensor, dense, displace, dagger, SparseOperator

function ps(ϕ::Float64,dim::FockBasis)
    x = -ϕ*im*number(dim)
    op = exp(dense(x))
    return op
end

function sms(g::Complex,dim::FockBasis)
    x = .5*(conj(g)*(destroy(dim)^2)-g*(create(dim)^2))
    op = exp(dense(x))
    return op
end

sms_re(g::Float64,dim::FockBasis) = sms(Complex(g),dim)
sms_im(g::Float64,dim::FockBasis) = sms(g*im,dim)

function disp(α::Complex,dim::FockBasis)
    op = displace(dim,α)
    return op
end

disp_re(α::Float64,dim::FockBasis) = disp(Complex(α),dim)
disp_im(α::Float64,dim::FockBasis) = disp(α*im,dim)

function tms(θ::Float64,dim::FockBasis)
    x = tensor(create(dim),create(dim))-tensor(destroy(dim),destroy(dim))
    x *= θ
    op = exp(dense(x))
    return op
end

function bs(θ::Float64,dim::FockBasis)
    x = tensor(create(dim),destroy(dim))+tensor(destroy(dim),create(dim))
    x *= im*θ
    op = exp(dense(x))
    return op
end

function swap(dim::FockBasis)
	op = sparse(bs(π/2,dim))
	return op
end

function apply(od::OpticalDevice,state::FockState)
	idd = identityoperator(state.dim)
	if typeof(od.mode) == Int
		mode1,mode2 = [od.mode,od.mode]
		Δmode = 1
	elseif length(od.mode) == 2
		mode1,mode2 = od.mode
		Δmode = abs(mode2-mode1)
	else
		throw(">2 modes operations are not implemented")
	end
	gate = od.optdev(od.param,state.dim)
	if Δmode == 1
		up = mode1 == 1 ? () : (idd for i in 1:(mode1-1))
		down = mode2 == state.n_mode ? () : (idd for i in 1:(state.n_mode-mode2))
		gate_ = tensor(up...,gate,down...)
	else
		gate_ = tensor(gate,(idd for i in Δmode-1)...)
		swaps = Vector{Operator}()
		swap_ = swap(state.dim)
		dom,upm = od.mode[2]>od.mode[1] ? od.mode : sort(od.mode) 
		for m in 1:Δmode-1
			up_ = (idd for i in 1:Δmode-m)
			down_ = (idd for i in 1:(m-1))
			s = tensor(up_...,swap_,down_...)
			push!(swaps,s)
		end
		gate_ = *(swaps...,gate_,reverse(swaps)...)
	end
	if isa(state.ρ, Ket)
		return gate_*state.ρ
	else
		return gate_*state.ρ*dagger(gate_)
	end
end

function apply!(od::OpticalDevice,state::FockState)
	ρ_ = apply(od,state)
	state.ρ = ρ_
end
