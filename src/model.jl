"""
    Layer(thick::Float64, sld::ComplexF64, rough::Float64)

A description of a layer in a system. 

Parameters
----------
- `thick::Float64` : the layer thickness.
- `sld::ComplexF64` : layer scattering length density.
- `rough::Float64` : roughness with layer above.
"""
struct Layer
    thick::Float64
    sld::ComplexF64
    rough::Float64
end