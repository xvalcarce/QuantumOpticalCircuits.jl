import QuantumOptics: create, destroy, number, identityoperator, tensor, dense, displace, dagger, SparseOperator

# Single-mode optical devices

# Phase-shifter
function ps(ϕ::Float64,dim::FockBasis)
    x = ϕ*im*number(dim)
    op = exp(dense(x))
    return op
end

# Single-mode squeezer
function sms(z::Complex,dim::FockBasis)
    x = .5*(conj(z)*(destroy(dim)^2)-z*(create(dim)^2))
    op = exp(dense(x))
    return op
end

function sms(r::Float64,ϕ::Float64,dim::FockBasis)
	z = r*exp(ϕ*im)
	op = sms(z,dim)
end

sms_re(g::Float64,dim::FockBasis) = sms(Complex(g),dim)
sms_im(g::Float64,dim::FockBasis) = sms(g*im,dim)

# Displacement
function disp(α::Complex,dim::FockBasis)
    op = displace(dim,α)
    return op
end

disp_re(α::Float64,dim::FockBasis) = disp(Complex(α),dim)
disp_im(α::Float64,dim::FockBasis) = disp(α*im,dim)

# Two-mode optical devices

# Two-mode squeezer
function tms(r::Float64,ϕ::Float64,dim::FockBasis)
	z = r*exp(ϕ*im)
	x = z*tensor(create(dim),create(dim))-conj(z)*tensor(destroy(dim),destroy(dim))
    op = exp(dense(x))
    return op
end

# Beam-splitter
function bs(θ::Float64,dim::FockBasis)
    x = tensor(create(dim),destroy(dim))+tensor(destroy(dim),create(dim))
    x *= im*θ
    op = exp(dense(x))
    return op
end

# Swap gate
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
	gate = od.optdev(od.param...,state.dim)
	if Δmode == 1
		up = mode1 == 1 ? () : (idd for i in 1:(mode1-1))
		down = mode2 == state.n_mode ? () : (idd for i in 1:(state.n_mode-mode2))
		gate_ = tensor(up...,gate,down...)
	else
		dom,upm = od.mode[2]>od.mode[1] ? od.mode : sort(od.mode) 
		idd_up = (idd for i in 1:(od.mode[1]-1))
		idd_do = (idd for i in (od.mode[2]+1):state.n_mode)
		gate_ = tensor(gate,(idd for i in 1:(Δmode-1))...)
		swap_ = Vector{Operator}()
		swap_gate = swap(state.dim)
		for m in 0:Δmode-2
			idd_swap_up_ = (idd for i in 1:(Δmode-1-m))
			idd_swap_down_ = (idd for i in 1:m)
			s = tensor(idd_swap_up_...,swap_gate,idd_swap_down_...)
			push!(swap_,s)
		end
		gate_ = *(swap_...,gate_,reverse(swap_)...)
		gate_ = tensor(idd_up...,gate_,idd_do...)
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
