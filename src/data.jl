using Measurements, DelimitedFiles, Distributions, LinearAlgebra

"""
    Data(q::Array{Measurements.Measurement}, R::Array{Measurements.Measurement}, resolution::Array{Distributions.UnivariateDistribution}, filepath::String)

The `struct` for storing experimental datasets.

### Parameters
- `q::Array{Measurement}` : the q-vector values, including dq uncertainty.
- `R::Array{Measurement}` : the reflected intensity values, including uncertainty.
- `resolution::Array{Distributions.UnivariateDistribution}` : the resolution function at each q-vector.
- `filepath::String` : the file path for the file.
"""
struct Data
    q::Array{Measurement}
    R::Array{Measurement}
    resolution::Array{Distributions.UnivariateDistribution}
    name::String
end

"""
    read_data(filename::String; delim=nothing, dq::Float64=0.05, dR::Float64=0.1)

Read experimental data from a file and store in a `Jeff.Data` object.

### Parameters
- `filename::String` : the file path to be read in.
- `delim::AbstractChar` : the delimiting character in the input file. If the file is whitespace delimited, pass `nothing`. Defaults to `nothing`.
- `dq::Float64` : percentage q-vector uncertainty, if not present in file. Defaults to `5.`.
- `dR::Float64` : percentage reflected intensity uncertainty, if not present in file. Defaults to `10.`.

### Returns
- `::Jeff.Data` : a data object containing the relevant information.
"""
function read_data(filename::String; delim=nothing, dq::Float64=5., dR::Float64=10.)
    if delim === nothing
        x = readdlm(filename)
    else
        x = readdlm(filename, delim)
    end
    nrows, ncols = size(x)
    q = zeros(Measurement, nrows)
    R = zeros(Measurement, nrows)
    resolution = []
    dq /= 100
    dR /= 100
    if ncols == 2
        for i = 1:nrows
            q[i] = measurement(x[i, 1], x[i, 1] * dq)
            R[i] = measurement(x[i, 2], x[i, 2] * dR)
            push!(resolution, Normal(x[i, 1], x[i, 1] * dq))
        end
    elseif ncols == 3
        for i = 1:nrows
            q[i] = measurement(x[i, 1], x[i, 1] * dq)
            R[i] = measurement(x[i, 2], x[i, 3])
            push!(resolution, Normal(x[i, 1], x[i, 1] * dq))
        end
    elseif ncols == 4
        for i = 1:nrows
            q[i] = measurement(x[i, 1], x[i, 4])
            R[i] = measurement(x[i, 2], x[i, 3])
            push!(resolution, Normal(x[i, 1], x[i, 4]))
        end
    else
        for i = 1:nrows
            q[i] = measurement(x[i, 1], x[i, 4])
            R[i] = measurement(x[i, 2], x[i, 3])
            distribution = getfield(Distributions, Symbol(x[i, 5]))
            push!(resolution, distribution(x[i, 6:end]...))
        end
    end
    return Data(q, R, resolution, filename)
end

"""
    get_distribution(y::Array{Measurements.Measurement})

Return the probability distribution for the data.

### Parameters
- `y::Array{Measurements.Measurement}` : ordinate data (reflected intensity).

### Returns
- `::Distributions.MvNormal` : probabilistic description of the data.
"""
function get_distribution(x::Array{Measurement}, y::Array{Measurement})
    distributions = []
    for i = 1:size(x, 1)
        position = [Measurements.value(x[i]), Measurements.value(y[i])]
        uncertainty = [Measurements.uncertainty(x[i]), Measurements.uncertainty(y[i])]
        push!(distributions, MvNormal(position, Diagonal(uncertainty)))
    end
    return distributions
end

"""
    transform(y::Array{Measurements.Measurement}, f)

Perform some transformation on the ordinate axis.

### Parameters
- `y::Array{Measurements.Measurement}` : ordinate data (reflected intensity).
- `f` : transformation to perform. Defaults to natural logarithm.

### Returns
- `::Array{Measurements.Measurement}` : transformed ordinate.
"""
function transform(y::Array{Measurement}; f=x->log.(x))
    return f(y)
end
