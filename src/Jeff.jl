module Jeff
using DelimitedFiles

include("reflect.jl")

function example_layers()
    layers = Array(zeros(Float64, 12))
    layers[2] = 0.
    layers[3] = 0.
    layers[5] = 100
    layers[6] = 3.47
    layers[8] = 2
    layers[10] = 2.07
    layers[12] = 3

    l = transpose(reshape(layers, (4, floor(Int64, size(layers, 1) / 4))))
    return l
end

function example_data()
    data = readdlm(string(pwd(), Base.Filesystem.path_separator, "../test/theoretical.txt"))
    return data
end

end # module
