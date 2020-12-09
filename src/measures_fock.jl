import QuantumOptics: sparse, create, destroy, tensor, ⊗, identityoperator, expect, tr, ptrace, normalize

function p_noclick(state::FockState,mode::Int,η::Float64)
	idd = identityoperator(state.dim)
	noclick_op = sparse(sum([((-(1-η))^n*create(state.dim)^n*destroy(state.dim)^n)/factorial(n) for n in 0:state.dim.N]))
	up = mode == 1 ? () : (idd for i in 1:(mode-1))
	down = mode == state.n_mode ? () : (idd for i in 1:(state.n_mode-mode))
	noclick_op_ = tensor(up...,noclick_op,down...)
	p_nc = expect(noclick_op_,state.ρ)
	return real(p_nc)
end

function p_noclick!(state::FockState,mode::Int,η::Float64)
	idd = identityoperator(state.dim)
	noclick_op = sparse(sum([((-(1-η))^n*create(state.dim)^n*destroy(state.dim)^n)/factorial(n) for n in 0:state.dim.N]))
	up = mode == 1 ? () : (idd for i in 1:(mode-1))
	down = mode == state.n_mode ? () : (idd for i in 1:(state.n_mode-mode))
	noclick_op_ = tensor(up...,noclick_op,down...)
	ρ_ = isa(state.ρ,Ket) ? state.ρ ⊗ dagger(state.ρ) : state.ρ
	ρ_ = noclick_op_*ρ_*dagger(noclick_op_)
	p_nc = abs(real(tr(ρ_)))
	ρ_ = normalize(ρ_)
	state.ρ = ρ_
	ptrace!(state,mode)
	return p_nc
end

function p_click!(state::FockState,mode::Int,η::Float64)
	idd = identityoperator(state.dim)
	noclick_op = sparse(sum([((-(1-η))^n*create(state.dim)^n*destroy(state.dim)^n)/factorial(n) for n in 0:state.dim.N]))
	click_op = idd-noclick_op
	up = mode == 1 ? () : (idd for i in 1:(mode-1))
	down = mode == state.n_mode ? () : (idd for i in 1:(state.n_mode-mode))
	click_op_ = tensor(up...,click_op,down...)
	ρ_ = isa(state.ρ,Ket) ? state.ρ ⊗ dagger(state.ρ) : state.ρ
	ρ_ = click_op_*ρ_*dagger(click_op_)
	p_c = abs(real(tr(ρ_)))
	ρ_ = normalize(ρ_)
	state.ρ = ρ_
	ptrace!(state,mode)
	return p_c
end

herald_noclick!(state::FockState,mode::Int,η::Float64) = p_noclick!(state,mode,η)
herald_click!(state::FockState,mode::Int,η::Float64) = p_click!(state,mode,η)

function E(state::FockState,k::Vector{Int},η::Float64)
	throw("POVM not implemented")
end

function correlator(state::FockState,i::Int,j::Int,η::Float64)
	if i==j
		throw(ArgumentError("Modes i,j need to be different"))
	elseif i>j
		i,j = j,i
	end
	if state.n_mode > 2
		m = filter(e->e∉[i,j],collect(1:state.n_mode))
		ptrace!(state,m)
	end
	idd = identityoperator(state.dim)
	noclick_op = sparse(sum([((-(1-η))^n*create(state.dim)^n*destroy(state.dim)^n)/factorial(n) for n in 0:state.dim.N]))
	#p_nc Alice
	noclick_op_i = tensor(noclick_op,idd)
	p_nc_a = real(expect(noclick_op_i,state.ρ))
	#p_nc Bob
	noclick_op_j = tensor(idd,noclick_op)
	p_nc_b = real(expect(noclick_op_j,state.ρ))
	#p_nc Cond
	noclick_op_ij = tensor(noclick_op,noclick_op)
	p_nc_ab = real(expect(noclick_op_ij,state.ρ))
	corr = 1 - 2*(p_nc_a+p_nc_b) + 4*p_nc_ab
    return corr
end
