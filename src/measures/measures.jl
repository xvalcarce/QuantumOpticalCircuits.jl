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
	tol::Int
end

PhotonDetector(mode::Int;η=0.0,tol=8) = PhotonDetector(mode,η,tol)

struct Heralding <: Measure
	mode::Int
	η::Float64
	tol::Int
end

Heralding(mode::Int;η=0.0,tol=8) = Heralding(mode,η,tol)

struct Wigner <: Measure
	mode::Int
	α::Union{Vector{Float64},Vector{Vector{Float64}}}
end

Wigner(mode::Int,α::Complex) = Wigner(mode,[real(α),imag(α)])
Wigner(mode::Int,α::Vector{Complex}) = Wigner(mode,[[real(a),imag(a)] for a in α])

function (meas::PhotonDetector)(state::AbstractState)
	out = p_click(state,meas.mode,meas.η)
	r_out = round(out, digits=meas.tol)
	return r_out
end

function (meas::Heralding)(state::AbstractState)
	out = herald!(state,meas.mode,meas.η,meas.tol)
	@debug "Heralding with p = $out"
	return state
end

function (meas::Wigner)(state::AbstractState)
	out = wigner(state,meas.mode,meas.α)
	return out
end
