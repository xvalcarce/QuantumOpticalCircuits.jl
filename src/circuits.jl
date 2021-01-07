import Base: print,show

mutable struct Circuit
	modes::Int
	state::State
	gates::Vector{Gate}
	meas_output::Vector{Tuple}
	function Circuit(modes,state)
		gates = Vector{Gate}()
		meas_output = Vector{Tuple}()
		new(modes,state,gates,meas_output)
	end
end

GaussianCircuit(modes::Int) = Circuit(modes,vacuum_gaussian(modes))
PseudoGaussianCircuit(modes::Int) = Circuit(modes,vacuum_pseudogaussian(modes))
FockCircuit(modes::Int,cutoff::Int) = Circuit(modes,vacuum_fock(modes,cutoff))

#TODO: Constructor from gates

function add_gate!(circuit::Circuit,gate::Gate)
	push!(circuit.gates,gate)
end

function add_gate!(circuit::Circuit,gate::Vector{<:Gate})
	append!(circuit.gates,gate)
end

function compute(circuit::Circuit;verbose=false)
	meas_output = Vector{Tuple}()
	state = copy(circuit.state)
	for (depth,gate) in enumerate(circuit.gates)
		if typeof(gate) == Measure
			r = apply!(gate,state)
			mo = (depth,String(Symbol(gate)),r)
			if verbose
				println(mo)
			end
			push!(meas_output,mo)
		else
			apply!(gate,state)
		end
	end
	return state,meas_output
end

function compute!(circuit::Circuit;verbose=false)
	for (depth,gate) in enumerate(circuit.gates)
		if typeof(gate) == Measure
			r = apply!(gate,circuit.state)
			mo = (depth,String(Symbol(gate)),r)
			if verbose
				println(mo)
			end
			push!(circuit.meas_output,mo)
		else
			apply!(gate,circuit.state)
		end
	end
end

function print(cir::Circuit)
	height = 2cir.modes+1
	asciir = Vector{String}()
	for i in 1:height
		if i%2 == 1
			push!(asciir,"   ")
		else
			push!(asciir,string(Int(i/2))*": ")
		end
	end
	bs_ascii = ["┤╲ ╱├","│ ╳ │","┤╱ ╲├"]
	tms_ascii = ["┤ T ├","│ M │","┤ S ├"]
	m_idx = 1
	for g in cir.gates
		if typeof(g.mode) == Int
			act = collect(2g.mode-1:2g.mode+1)
			idle = filter(e->e∉act,collect(1:height))
			for i in idle
				if i%2 == 1
					asciir[i] *= "     "
				else
					asciir[i] *= "─────"
				end
			end
			if typeof(g) == Measure
				asciir[act[1]] *= "     "
				asciir[act[2]] *= "──◗~ "
				asciir[act[3]] *= "     "
				if length(cir.meas_output) != 0
					asciir[act[2]] *=  string(cir.meas_output[m_idx][3])
					m_idx += 1
				end
				asciir[act[2]] *= "▢"
			else
				sym = uppercase(string(Symbol(g.optdev))[1]) 
				asciir[act[1]] *= "┌───┐" 
				asciir[act[2]] *= "┤ "*sym*" ├"
				asciir[act[3]] *= "└───┘"
			end
		else
			down,up = sort(g.mode)
			act = collect(2down-1:2up+1)
			idle = filter(e->e∉act,collect(1:height))
			for i in idle
				if i%2 == 1
					asciir[i] *= "     "
				else
					asciir[i] *= "─────"
				end
			end
			asciir[act[1]] *= "┌───┐" 
			asciir[act[end]] *= "└───┘"
			if g.optdev ∈ [bs,swap]
				gate_ascii = bs_ascii
			elseif g.optdev == tms
				gate_ascii = tms_ascii
			else
				throw("Printing gate "*string(Symbol(g.optdev))*" is not implemented")
			end
			if up-down == 1
				for (idx,a) in enumerate(act[2]:act[4])
					asciir[a] *= gate_ascii[idx]
				end
			else
				mid = act[Int((length(act)+1)/2)]
				asciir[act[2]] *= gate_ascii[1]
				for a in act[2]+1:mid-1
					asciir[a] *= "┊   ┊"
				end
				asciir[act[mid]] *= gate_ascii[2]
				for a in mid+1:length(asciir)-2
					asciir[a] *= "┊   ┊"
				end
				asciir[act[end-1]] *= gate_ascii[3]
			end
		end
	end
	la = maximum(map(length,asciir))
	for (idx,a) in enumerate(asciir)
		if findfirst("▢",a) != nothing
			a_, _ = split(a,"▢")
			a_ *= " "^(la-length(a)+1)
			asciir[idx] = a_
		end
	end
	la = maximum(map(length,asciir))
	for (idx,a) in enumerate(asciir)
		asciir[idx] *= " "^(la-length(a))
	end
	return asciir
end

Base.show(io::IO,cir::Circuit) = display(print(cir))
