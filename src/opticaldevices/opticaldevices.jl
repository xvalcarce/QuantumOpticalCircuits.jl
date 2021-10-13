include("./operators/operators.jl")

abstract type AbstractOpticalDevice <: Gate end

struct SingleModeOpticalDevice <: AbstractOpticalDevice
	mode::Int
	op::SingleModeOperator
end

struct TwoModeOpticalDevice <: AbstractOpticalDevice 
	modes::Vector{Int}
	op::TwoModeOperator
end

# Single-mode optical devices

# Phase shifter
PhaseShifter(mode::Int,param::Float64) = SingleModeOpticalDevice(mode,PhaseShifter(param))

# Single-mode squeezer
# z = r e^{iϕ}
SMS(mode::Int,z::Complex) = SingleModeOpticalDevice(mode,SMS(z))
SMS(mode::Int,r::Float64,ϕ=0.0) = SingleOpticalDevice(mode,SMS(r,ϕ))

# Displacement
# α = r e^{iϕ}
D(mode::Int,α::Complex) = VariableOpticalDevice(mode,α,disp)
D(mode::Int,r::Float64,ϕ=0.0) = VariableOpticalDevice(mode,r*exp(ϕ*im),disp)
D_Re(mode::Int,param::Float64) = VariableOpticalDevice(mode,param,disp_re)
D_Im(mode::Int,param::Float64) = VariableOpticalDevice(mode,param,disp_im)

# Two-mode optical devices

# Two-mode sequeezer
# z = r e^{iϕ}
TMS(mode::Vector{Int},r::Float64,ϕ=0.0) = VariableOpticalDevice(mode,[r,ϕ],tms)

# Beam-splitter
# default to 50-50 BS
BS(mode::Vector{Int},θ=π/4) = VariableOpticalDevice(mode,θ,bs)

# Swap gate
SWAP(mode::Vector{Int}) = ConstantOpticalDevice(mode,swap)

# Apply

function apply(od::AbstractOpticalDevice,state::State)
	state_ = copy(state)
	apply!(od,state_)
	return state_
end


