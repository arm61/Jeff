# reflect.jl

These functions are focused on the calculation of a reflectometry profile from an abstract set of layers.

```@docs
Jeff.constant_smearing(q::Array{Float64, 1}, w::Array{Any, 2} resolution::Any, scale::Any, bkg::Any)
Jeff.abeles(q::Array{Float64, 1}, layers::Array{Float64, 2})
Jeff.same_convolution(a::Array{Float64, 1}, b::Array{Float64, 1})
Jeff.gauss(x::Array{Float64, 1}, s::Float64)
```
