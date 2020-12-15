struct OpticalDevice <: Gate
	mode::Union{Int,Array{Int,1}}
	param::Union{Float64,Complex}
	optdev::Function
end

include("opticaldevices_gaussian.jl")
include("opticaldevices_fock.jl")

PS(mode::Int,param::Float64) = OpticalDevice(mode,param,ps)

SMS(mode::Int,param::Complex) = OpticalDevice(mode,param,sms)
SMS(mode::Int,param::Float64) = OpticalDevice(mode,param,sms_re)
SMS_Re(mode::Int,param::Float64) = OpticalDevice(mode,param,sms_re)
SMS_Im(mode::Int,param::Float64) = OpticalDevice(mode,param,sms_im)

D(mode::Int,param::Complex) = OpticalDevice(mode,param,disp)
D(mode::Int,param::Float64) = OpticalDevice(mode,param,disp_re)
D_Re(mode::Int,param::Float64) = OpticalDevice(mode,param,disp_re)
D_Im(mode::Int,param::Float64) = OpticalDevice(mode,param,disp_im)

TMS(mode::Array{Int,1},param::Float64) = OpticalDevice(mode,param,tms)

BS(mode::Array{Int,1},param::Float64) = OpticalDevice(mode,param,bs)
SWAP(mode::Array{Int,1}) = OpticalDevice(mode,swap)
