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
FockCircuit(modes::Int,cutoff::Int) = Circuit(modes,vacuum_fock(modes,cutoff))

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
	for (depth,gate) in circuit.gates
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
