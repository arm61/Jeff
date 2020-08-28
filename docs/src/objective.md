# objective.jl

The code here is focused on optimisation and sampling

```@docs
Jeff.variables_and_bounds(model::Jeff.Model)
Jeff.populate_model(variables::Array{Float64, 1}, model::Jeff.Model)
Jeff.loglikelihood(distributions, model::Array{Float64, 2})
Jeff.negativeloglikelihood(distributions, model::Array{Float64, 2})
```
