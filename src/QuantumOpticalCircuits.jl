module QuantumOpticalCircuits

export State,
	GaussianState, vacuum_gaussian,
	PseudoGaussianState, vacuum_pseudogaussian,
	FockState, vacuum_fock,
	OpticalDevice, apply, apply!,
	PS, SMS, SMS_Re, SMS_Im, D, D_Re, D_Im, TMS, BS,
	Measure, correlator, 
	PhotonDetector, Heralding,
	Circuit, FockCircuit, GaussianCircuit, PseudoGaussianCircuit,
	add_gate!, compute, compute!, print

include("states/states.jl")

abstract type Gate end

include("opticaldevices/opticaldevices.jl")
include("measures/measures.jl")

include("circuits.jl")

end
