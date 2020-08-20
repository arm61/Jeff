using Measurements, DelimitedFiles

struct Data
    q::Array{Measurement}
    R::Array{Measurement}
end

function read_data(filename::String, delim=nothing, dq=0.05, dR=0.1)
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
    return Data(q, R)
end