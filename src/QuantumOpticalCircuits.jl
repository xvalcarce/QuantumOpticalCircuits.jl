module QuantumOpticalCircuits

export State,
	GaussianState, vacuum_gaussian,
	FockState, vacuum_fock,
	OpticalDevice, apply, apply!,
	PS, SMS, SMS_Re, SMS_Im, D, D_Re, D_Im, TMS, BS,
	Measure, correlator, 
	PhotonDetection, Heralding,
	Circuit, FockCircuit, GaussianCircuit, 
	add_gate!, compute, compute!

include("states.jl")
include("states_gaussian.jl")
include("states_fock.jl")

abstract type Gate end
include("opticaldevices.jl")
include("measures.jl")

include("circuits.jl")

end
