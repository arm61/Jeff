using Interpolations

const TINY = 1e-30
const _FWHM = 2 * sqrt(2 * log(2.0))

function abeles(q_vals::Array{Float64, 1}, layers::Array{Float64, 2}, scale=1., bkg=0.)
    nlayers = size(layers, 1) - 2
    q = hcat(fill(q_vals, size(layers, 1))...)
    npnts = size(q, 1)
    
    mi00 = ones(ComplexF64, npnts, nlayers + 1)
    
    sld = zeros(ComplexF64, 1, nlayers + 2)
    
    sld[1, 2:end] += ((layers[2:end, 2] .- layers[1, 2]) + 1im * (abs.(layers[2:end, 3]) .+ TINY)) * 1e-6
            
    kn = sqrt.(q .^ 2 ./ 4 .- 4 .* pi .* sld)::Array{ComplexF64, 2}
        
    rj = kn[:, 1:end-1] .- kn[:, 2:end]
    rj ./= (kn[:, 1:end-1] + kn[:, 2:end])
    rj .*= exp.(-2 .* kn[:, 1:end-1] .* kn[:, 2:end] .* transpose(layers[2:end, 4]) .^ 2)
    
    if nlayers > 0
        mi00[:, 2:end] = exp.(kn[:, 2:end-1] * 1im * abs.(layers[2:end-1, 1]))
    end
    mi11 = 1. ./ mi00
    mi10 = rj .* mi00
    mi01 = rj .* mi11    
    
    mrtot00 = mi00[:, 1]
    mrtot01 = mi01[:, 1]
    mrtot10 = mi10[:, 1]
    mrtot11 = mi11[:, 1]
    
    for idx in range(2, nlayers + 1, step=1)
        p0 = mrtot00 .* mi00[:, idx] .+ mrtot10 .* mi01[:, idx]
        p1 = mrtot00 .* mi10[:, idx] .+ mrtot10 .* mi11[:, idx]
        mrtot00 = p0
        mrtot10 = p1

        p0 = mrtot01 .* mi00[:, idx] .+ mrtot11 .* mi01[:, idx]
        p1 = mrtot01 .* mi10[:, idx] .+ mrtot11 .* mi11[:, idx]

        mrtot01 = p0
        mrtot11 = p1
    end
            
    r = (mrtot01 ./ mrtot00)
    
    reflectivity = r .* conj(r)
    reflectivity *= scale
    reflectivity .+= bkg
    reflectivity = real(reflectivity)
    return reflectivity
end

function same_convolution(a, b)
    m = length(a)
    n = length(b)
    
    padleft = ceil(Int32, n/2) - 1
    padright = floor(Int32, n/2)
    
    a = append!(zeros(padleft), a)
    a = append!(a, zeros(padright))
    
    output = zeros(m)
    for i in 1:m
        for j in 1:n
            output[i] += a[i+j-1] * b[n-j+1]
        end
    end
    return output
end

function gauss(x, s)
    g = 1. / s / sqrt(2 * pi) * exp.(-0.5 * x .^ 2 / s / s)
    return g
end

function constant_smearing(q::Array{Float64, 1}, w::Array{Float64, 2}, resolution)
    if resolution < 0.5
        return abeles(q, w)
    end

    resolution /= 100
    gaussnum = 51
    gaussgpoint = (gaussnum - 1) / 2
    

    lowq = minimum(q[:])
    highq = maximum(q[:])
    if lowq <= 0.
        lowq = 1e-6
    end

    start = log10(lowq) - 6 * resolution / _FWHM
    finish = log10(highq * (1 + 6 * resolution / _FWHM))
    interpnum = round(abs(1 * (abs(start - finish)) / (1.7 * resolution / _FWHM / gaussgpoint)))
    xtemp = range(start, finish, length=Int(interpnum))
    xlin = 10. .^ xtemp
    
    gauss_x = range(-1.7 * resolution, 1.7 * resolution, length=gaussnum)
    gauss_y = gauss(gauss_x, resolution / _FWHM)

    rvals = abeles(xlin, w)
    smeared_rvals = same_convolution(rvals, gauss_y)
    smeared_rvals *= gauss_x[2] - gauss_x[1]        
 
    itp = LinearInterpolation(xlin, smeared_rvals)
    smeared_output = itp(q)
    return smeared_output
end