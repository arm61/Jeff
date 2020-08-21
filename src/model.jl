"""
    Parameter(value::Float64, vary::Bool, prior::UnivariateDistribution)

A parameter that can be optimised in the analysis procress.

Parameters
----------
- `value::Float64` : the value for the given parameter.
- `vary::Bool` : should the parameter be varied in optimisation. 
- `prior::Distrbutions.UnivariateDistribution` : the priopr probability distribution for the parameter. 
"""
struct Parameter
    value::Float64
    vary::Bool
    prior
end

"""
    Layer(thick::Jeff.Parameter, sld::Jeff.Parameter, isld::Jeff.Parameter, rough::Jeff.Parameter)

A description of a layer in a system. 

Parameters
----------
- `thick::Jeff.Parameter` : the layer thickness.
- `sld::Jeff.Parameter` : layer real scattering length density.
- `isld::Jeff.Parameter` : layer imaginary scattering length density.
- `rough::Jeff.Parameter` : roughness with layer above.
"""
struct Layer
    thick::Jeff.Parameter
    sld::Jeff.Parameter
    isld::Jeff.Parameter
    rough::Jeff.Parameter
end