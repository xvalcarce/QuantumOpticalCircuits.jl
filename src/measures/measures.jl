export PhotonDetector, Heralding, Wigner

include("measures_gaussian.jl")
include("measures_fock.jl")

function p_click(state::AbstractState,mode::Int,η::Float64)
	p_nc = p_noclick(state,mode,η)
	p_c = 1-p_nc
	return p_c
end

abstract type Measure end

struct PhotonDetector <: Measure
	mode::Int
	η::Float64
	tol::Float64
end

struct Heralding <: Measure
	mode::Int
	η::Float64
	tol::Float64
end

struct HeraldingNoClick <: Measure
	mode::Int
	η::Float64
	tol::Float64
end

struct Wigner <: Measure
	α::Union{Vector{Float64},Vector{Vector{Float64}}}
end

function (meas::PhotonDetector)(state::AbstractState)
	out = p_click(state,meas.mode,meas.η)
    r_out = out > meas.tol ? out : 0.0
	return r_out
end

function (meas::Heralding)(state::AbstractState)
	out = herald_click!(state,meas.mode,meas.η,meas.tol)
	return state
end

function (meas::HeraldingNoClick)(state::AbstractState)
	out = herald_noclick!(state,meas.mode,meas.η)
    @assert out > meas.tol "NoClick probability below tolerance ($meas.tol): $out"
	return state
end

function (meas::Wigner)(state::AbstractState)
	out = wigner(state,meas.α)
	return out
end

PhotonDetector(mode::Int;η=1.0,tol=1e-12) = PhotonDetector(mode,η,tol)
Heralding(mode::Int;η=1.0,tol=1e-12,onclick=true) = onclick ? Heralding(mode,η,tol) : HeraldingNoClick(mode,η,tol)
