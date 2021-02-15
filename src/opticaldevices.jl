struct OpticalDevice <: Gate
	mode::Union{Int,Vector{Int}}
	param::Union{Float64,Complex,Vector{Float64}}
	optdev::Function
end

include("opticaldevices_gaussian.jl")
include("opticaldevices_fock.jl")

# Single-mode optical devices

# Phase shifter
PS(mode::Int,param::Float64) = OpticalDevice(mode,param,ps)

# Single-mode squeezer
# z = r e^{iϕ}
SMS(mode::Int,z::Complex) = OpticalDevice(mode,z,sms)
SMS(mode::Int,r::Float64,ϕ=0.0) = OpticalDevice(mode,[r,ϕ],sms)
SMS_Re(mode::Int,param::Float64) = OpticalDevice(mode,param,sms_re)
SMS_Im(mode::Int,param::Float64) = OpticalDevice(mode,param,sms_im)

# Displacement
# α = r e^{iϕ}
D(mode::Int,α::Complex) = OpticalDevice(mode,α,disp)
D(mode::Int,r::Float64,ϕ=0.0) = OpticalDevice(mode,r*exp(ϕ*im),disp)
D_Re(mode::Int,param::Float64) = OpticalDevice(mode,param,disp_re)
D_Im(mode::Int,param::Float64) = OpticalDevice(mode,param,disp_im)

# Two-mode optical devices

# Two-mode sequeezer
# z = r e^{iϕ}
TMS(mode::Vector{Int},r::Float64,ϕ=0.0) = OpticalDevice(mode,[r,ϕ],tms)

# Beam-splitter
# default to 50-50 BS
BS(mode::Vector{Int},θ=π/4) = OpticalDevice(mode,θ,bs)

# Swap gate
SWAP(mode::Vector{Int}) = OpticalDevice(mode,swap)
