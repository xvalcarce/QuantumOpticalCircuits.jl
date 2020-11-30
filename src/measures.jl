include("measures_gaussian.jl")
include("measures_fock.jl")

function p_click(state::State,mode::Int,η::Float64)
	p_nc = p_noclick(state,mode,η)
	p_c = 1-p_nc
	return p_c
end

function p_click!(state::State,mode::Int,η::Float64)
	p_nc = p_noclick!(state,mode,η)
	p_c = 1-p_nc
	return p_c
end

struct Measure <: Gate
	mode::Int
	measdev::Function
	η::Float64
end

function apply!(meas::Measure,state::State)
	res = meas.measdev(state,meas.mode,meas.η)
	return res
end

PhotonDetection(mode::Int;η=0.0) = Measure(mode,p_click,η)
Heralding(mode::Int;η=0.0) = Measure(mode,p_click!,η)
