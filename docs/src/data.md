# data.jl

Functions to import experimental datasets and store them as [`Jeff.Data`](@ref) objects. 

```@docs
Jeff.read_data(filename::String, delim=nothing, dq::Float64=0.05, dR::Float64=0.1)
Jeff.Data
```