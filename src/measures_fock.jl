import QuantumOptics: sparse, create, destroy, tensor, ⊗, identityoperator, expect, tr, ptrace

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
	ρ_ /= tr(ρ_)
	state.ρ = ptrace(ρ_,mode)
	state.n_mode -= 1
	return p_nc
end

function correlator(state::FockState,i::Int,j::Int,η::Float64)
	if i==j
		throw(ArgumentError("Modes i,j need to be different"))
	elseif i>j
		i,j = j,i
	end
	idd = identityoperator(state.dim)
	noclick_op = sparse(sum([((-(1-η))^n*create(state.dim)^n*destroy(state.dim)^n)/factorial(n) for n in 0:state.dim.N]))
	#p_nc Alice
	up_i = i == 1 ? () : (idd for i in 1:(i-1))
	down_i = (idd for i in 1:(state.n_mode-i))
	noclick_op_i = tensor(up_i...,noclick_op,down_i...)
	p_nc_a = real(expect(noclick_op_i,state.ρ))
	#p_nc Bob
	up_j = (idd for i in 1:(j-1))
	down_j = j == state.n_mode ? () : (idd for i in 1:(state.n_mode-j))
	noclick_op_j = tensor(up_j...,noclick_op,down_j...)
	p_nc_b = real(expect(noclick_op_j,state.ρ))
	#p_nc Cond
	i,j = i<j ? (i,j) : (j,i)
	mid_ij = j-i == 1 ? () : (idd for i in 1:((j-i)-1))
	noclick_op_ij = tensor(up_i...,noclick_op,mid_ij...,noclick_op,down_j...)
	p_nc_ab = real(expect(noclick_op_ij,state.ρ))
	corr = 1 - 2*(p_nc_a+p_nc_b) + 4*p_nc_ab
    return corr
end
