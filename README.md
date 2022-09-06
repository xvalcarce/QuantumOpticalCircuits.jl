# QuantumOpticalCircuits.jl

**QuantumOpticalCircuits.jl** is a numerical framework for the simulation of parametric quantum optical circuits.
Written in [JULIA](https://julialang.org), this package aim for efficient and accurate quantum optics simulations, while keeping simplicity of use at its core.

## 🎁 Installation

This package is still in developpement and is thus not yet available in Julia's General Regisitries.
Installing is done using git

```julia
pkg> add https://github.com/xvalcarce/QuantumOpticalCircuits.jl
```

## 👀 Overview

`QuantumOpticalCircuits.jl` simulate quanutm optical circuits from 3 components: `States`, `OpticalDevices`, and `Measures`.  

A `n`-mode state can be a `FockState(n,d)`, a `GaussianState(n)` or a `PseudoGaussianState(n)`. The first type uses the Fock representation where `d` is the maximum number of photons considered. The last two type uses a Gaussian representation. `PseudoGaussianState` allows for non-Gaussian operation such as heralding operations (this is done by allowing the state to be a linear combination of Gaussian states).

Optical devices are `Operators` applied on specific optical modes. `Operators` implemented are phase-shifters (`PS`), displacements (`D`), single-mode squeezers (`SMS`), two-mode squeezers (`TMS`) and beam-splitters (`BS`).

Measurement are embbed in a type that can be applied on any state. We provide non-photon number resolving photon detector (`PhotonDetector`) and heralding operation (`Heradling`). Homodyne measurement are WIP.

Simple example:
```julia
julia> using QuantumOpticalCircuits
# 2-mode vacuum state using Gaussian represenation
julia> state = GaussianState(2)
# Let's apply a squeezing operation on mode 2, followed by a 50:50 beam-splitter between mode 1 and 2
julia> state = state |> SMS(0.42)(2) |> BS(π/4)(1,2)
# Finally, let's measure the probability to get a click in the 1st mode with a photon detector with 80% efficiency
julia> state |> PhotonDetector(1,η=0.8)
0.05495903420401427
# Or in one-line
julia> GaussianState(2) |> SMS(0.42)(2) |> BS(π/4)(1,2) |> PhotonDetector(1,η=0.8)
0.05495903420401427
# Or using the Fock representation with a cut-off at 5 photons
julia> FockState(2,5) |> SMS(0.42)(2) |> BS(π/4)(1,2) |> PhotonDetector(1,η=0.8)
0.054904814588312534
```

More example are available in the [examples](https://github.com/xvalcarce/QuantumOpticalCircuits.jl/tree/master/examples) folder.


## 🙏 Citation & Support

If you find `QuantumOpticalCircuits.jl` helpful, or if you simply want to support the project, we would appreciate if you starred the project!  
Contribution in any form (submitting issue, enhancement suggestion, pull requests...) are more welcome.  
To cite this project, please use the informations available in the [`CITATION.bib`](https://github.com/xvalcarce/QuantumOpticalCircuits.jl/blob/master/CITATION.bib) file.
