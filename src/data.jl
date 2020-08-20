using Measurements, DelimitedFiles

"""
    Data(q::Array{Measurements.Measurement}, R::Array{Measurements.Measurement})

The `struct` for storing experimental datasets. 

Parameters
----------
- `q::Array{Measurement}` : the q-vector values, including dq uncertainty.
- `R::Array{Measurement}` : the reflected intensity values, including uncertainty.
- `filepath::String` : the file path for the file. 
"""
struct Data
    q::Array{Measurement}
    R::Array{Measurement}
    filepath::String
end

"""
    read_data(filename::String, delim=nothing, dq::Float64=0.05, dR::Float64=0.1)

Read experimental data from a file and store in a `Jeff.Data` object. 

Parameters
----------
- `filename::String` : the file path to be read in. 
- `delim::AbstractChar` : the delimiting character in the input file. If the file is whitespace delimited, pass `nothing`. Defaults to `nothing`.
- `dq::Float64` : percentage q-vector uncertainty, if not present in file. Defaults to `0.05`.
- `dR::Float64` : percentage reflected intensity uncertainty, if not present in file. Defaults to `0.1`.

Returns
-------
- `::Jeff.Data` : a data object containing the relevant information.
"""
function read_data(filename::String, delim=nothing, dq::Float64=0.05, dR::Float64=0.1)
    if delim == nothing
        x = readdlm(filename)
    else
        x = readdlm(filename, delim)
    end
    nrows, ncols = size(x)
    q = zeros(Measurement, nrows)
    R = zeros(Measurement, nrows)
    if ncols == 2
        for i in range(1, nrows, step=1)
            q[i] = measurement(x[i, 1], x[i, 1] * dq)
            R[i] = measurement(x[i, 2], x[i, 2] * dR)
        end
    elseif ncols == 3
        for i in range(1, nrows, step=1)
            q[i] = measurement(x[i, 1], x[i, 1] * dq)
            R[i] = measurement(x[i, 2], x[i, 3])
        end
    elseif ncols == 4
        for i in range(1, nrows, step=1)
            q[i] = measurement(x[i, 1], x[i, 4])
            R[i] = measurement(x[i, 2], x[i, 3])
        end
    else
        throw(ArgumentError("The file that you are opening does not have 2, 3, or 4 columns"))
    end
    return Data(q, R, filename)
end