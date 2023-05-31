using Pkg

try
Pkg.activate("dimensionality_pkg")
catch
end
using Dates
using FileIO
using JLD2
using Statistics
using StatsBase
using MultivariateStats
using Random
using CSV
using DataFrames
using FFTW
using LinearAlgebra
using ManifoldLearning
using FourierAnalysis
using ProgressMeter
using ProfileSVG
using LsqFit
using FileIO
using ImageIO
using NearestNeighbors
using SharedArrays

# This function initializes logging by creating a directory for the logs and returning the path to the directory.

function init_logging()
    # Get the current week number based on the current unix time.
    weekno = week(unix2datetime(time()))

    # Create a directory string with the format "KW_2_<week number with leading zero>/"
    datestring = string("KW_2_", lpad(weekno, 2, "0"), "/")

    # Set the parent directory where logs will be stored.
    workdir = "/net/home/lschulz/logs/"

    # Append the datestring to the parent directory to create the full directory path.
    dir = workdir * datestring

    # Create the directory if it doesn't exist already.
    if isdir(dir) == false
        mkdir(dir)
    end

    # Return the directory path.
    return dir
end



# This function standardizes the input data by subtracting the mean and dividing by the standard deviation.

function centralizer(data::Vector{Float32})
    # Calculate the mean and standard deviation of the input data.
    m = mean(data)
    cv = std(data)

    # Subtract the mean from each element of the input data and divide by the standard deviation.
    return (data .- m) ./ cv
end


# This function embeds a lagged version of the input data into a matrix.

function embed_lag(data::Vector{Float32}, W::Int64)
    # Initialize an empty vector Y to store the embedded data.
    Y = []

    # Loop through the data and embed a lagged version of the data into Y.
    for i = 1:length(data) - W + 1
        Y = append!(Y, data[i:i + W - 1])
    end

    # Convert the embedded data to a matrix with W rows and (length(data) - W + 1) columns, and return it.
    return reshape(float.(Y), W, length(data) - W + 1)::Matrix{Float32}
end


# This function reconstructs a time series from a given set of principal components and eigenvectors.

function reconstructor(A, rho, N, M)
    # Calculate the length of the reconstructed time series.
    P = N - M + 1
    
    # Define a function R to compute the reconstructed values for a given time t, based on the PC projection A and eigenvector rho.
    R(t, M_t, L_t, U_t) = M_t * sum([A[t - j + 1] * rho[j] for j in L_t:U_t])
    
    # Define a function choicer to choose the parameters for R based on the value of t.
    function choicer(t)
        if 1 <= t <= M - 1
            return 1 / t, 1, t
        elseif M <= t <= P
            return 1 / M, 1, M
        elseif P + 1 <= t <= N
            return 1 / (N - t + 1), t - N + M, M
        end
    end
    
    # Loop through each time t and compute the reconstructed value using the R function and the choicer function.
    return [R(t, choicer(t)...) for t in 1:N]
end


"""
the calculation
"""

mutable struct matrixholder
    #requirements for function that puts in data
    N           ::Int64    # Number of data points
    W           ::Int64    # Embedding window length
    k           ::Int64    # Number of modes to keep
    signal      ::Vector{Float32}    # Original signal of length N
    emb         ::Matrix{Float32}    # Embedded signal of size (N-W+1) x W

    #requirements for function that finds out epsilon
    eps         ::Int64    # Epsilon value

    #requirements for dimensionality reduction
    lambda      ::Vector{Float32}    # Eigenvalues of EOF
    EOF         ::Matrix{Float32}    # Empirical Orthogonal Functions of size W x k
    PC         ::Matrix{Float32}    # Principal Components of size (N-W+1) x k
    RC         ::Matrix{Float32}    # Reconstructed signal of size N x k

    #the init function
    function matrixholder(N,W,k)
        # Initialize variables with undefined values
        signal      = Vector{Float32}(undef,N)
        emb         = Matrix{Float32}(undef,N-W+1,W)
        eps         = 0
        EOF         = Matrix{Float32}(undef,W,k)
        lambda      = Vector{Float32}(undef,k)
        PC         = Matrix{Float32}(undef,N-W+1,k)
        RC         = Matrix{Float32}(undef,N,k)
        return new(
            N,
            W,
            k,
            signal,
            emb,
            eps,
            lambda,
            EOF,
            PC,
            RC
        )
    end
end


function put_data_gSSA(d::matrixholder,signal::Vector{Float32},sampleyear::Int64)
    d.signal = centralizer(deseasonalize_gSSA(centralizer(signal),sampleyear::Int64))
    d.emb = centralized_embed_lag(d.signal,d.W)'
end

function put_data_lSSA(d::matrixholder,signal::Vector{Float32},sampleyear::Int64)
    d.signal = centralizer(deseasonalize_lSSA(centralizer(signal),sampleyear::Int64))
    d.emb = centralized_embed_lag(d.signal,d.W)'
end

function put_data_window(d::matrixholder,signal::Vector{Float32},sampleyear::Int64)
    windowlength=Int(floor(sampleyear/6))
    overlap=Int(floor(sampleyear/12))
    d.signal = centralizer(moving_window_deseason(centralizer(signal),windowlength,overlap))
    d.emb = centralized_embed_lag(d.signal,d.W)'
end

function put_data_raw(d::matrixholder,signal::Vector{Float32})
    d.signal = centralizer(signal)
    d.emb = centralized_embed_lag(d.signal,d.W)'

end

function centralized_embed_lag(data::Vector{Float32},W::Int64)
    #centralize
    function centralizer(data)
        m = mean(data)
        cv = std(data)
        return  (data.-m)./cv
    end
    Y = []
    for i=1:length(data)-W+1
        Y = append!(Y,centralizer(data[i:i+W-1]))
    end
    return reshape(float.(Y),W,length(data)-W+1)::Matrix{Float32}
end

"""
#the epsilon findout
"""
#required struct to iterate with different epsilon
mutable struct epsilon
    data_samples :: Matrix{Float32}
    eps         :: Vector{Float32}
    L           :: Vector{Float32}
    Weight      :: Matrix{Float32}
    function epsilon(W,stepnumber,sampling_size)
        min_eps = Float32(10^0)
        max_eps = Float32(10^7)

        eps             = 10 .^ range(log10(min_eps),log10(max_eps),length=stepnumber)
        data_samples    = Matrix{Float32}(undef,W,sampling_size)
        L               = Vector{Float32}(undef,length(eps))
        Weight          = Matrix{Float32}(undef,sampling_size,sampling_size)

        return new(data_samples,eps,L,Weight)
    end
end

# iteration function that returns the median
function fit_epsilon(data::Matrix{Float32},stepnumber,sampling_size,it_number,W)
    P = size(data)[2]
    object = epsilon(W,stepnumber,sampling_size)
    fit_eps_L = Vector{Float32}(undef,it_number)
    sample = rand(1:P,sampling_size,it_number)
    for t in 1:it_number
        object.data_samples = data[:,sample[:,t]]

        for (i,eps) in enumerate(object.eps)
            for i in 1:sampling_size, j in 1:sampling_size
                object.Weight[i,j] = exp(- norm(object.data_samples[:,i] - object.data_samples[:,j])^2 / eps)
            end
            object.L[i] = sum(object.Weight)

        end

        p0 = ones(3)
        model(eps,p) = p[1] .* atan.(eps .- p[2]) .+ p[3]
        p_midpoint = coef(curve_fit(model, log10.(object.eps), log10.(object.L), p0))
        a = p_midpoint[2]
        b = median(log10.(object.L[1:8]))
        one_over_e = model(b+(a-b)/exp(1),p_midpoint)
        fit_eps_L[t] = 10^(one_over_e)
    end
    return Int64(floor(median(fit_eps_L)))
end

#function that works on object
function put_epsilon(d::matrixholder)
    stepnumber      = 32
    sampling_size   = 64
    it_number       = 4

    d.eps = fit_epsilon(d.emb,stepnumber,sampling_size,it_number,d.W)
end
"""
#functions that fit in the model
"""
#little struct to hold all the variable types for the fit?
struct diffholder
    f::DiffMap{Float32}
    function diffholder(data::Matrix{Float32},k,t,alpha,eps)
        return new(fit(DiffMap,data,maxoutdim=k,t=1, α=alpha, ɛ=eps))
    end
end

#type-stable extraction
function gettheproj_diff(diff::diffholder)
    return Matrix(diff.f.proj')::Matrix{Float32},diff.f.λ::Vector{Float32}
end

#function that works on the object                  # FIXED t=1 alpha = 1
function put_EOF_diff(d::matrixholder)

    t = 1
    alpha = 1.0
    eps = d.eps

    diff = diffholder(d.emb,d.k,t,alpha,eps)
    d.EOF,d.lambda = gettheproj_diff(diff)

end

struct ssaholder
    f::PCA{Float32}
    function ssaholder(data::Matrix{Float32},k)
        return new(fit(PCA,data',maxoutdim=k,method=:auto,pratio=1))
    end
end

#type-stable extraction
function gettheproj_ssa(ssa::ssaholder)
    return Matrix(ssa.f.proj)::Matrix{Float32},(ssa.f.prinvars ./ ssa.f.prinvars[1])::Vector{Float32}
end

#function that works on the object
function put_EOF_ssa(d::matrixholder)
    d.EOF = zeros(Float32,d.W,d.k)
    d.lambda = zeros(Float32,d.k)

    ssa = ssaholder(d.emb,d.k)
    EOF,lambda = gettheproj_ssa(ssa)

    d.EOF[:,1:size(EOF)[2]] = EOF
    d.lambda[1:size(lambda)[1]] = lambda

end

"""
#building PC,RC   
"""                  
function calculate(d::matrixholder)
    d.PC = d.emb * d.EOF # hcat([pc!(d.signal,d.EOF[:,i],d.N-d.W-1,d.W) for i in 1:d.k]...) # 
    d.RC = hcat([reconstructor(d.PC[:,i],d.EOF[:,i],d.N,d.W) for i in 1:d.k]...)
    d.lambda = Float32.(diag(1/d.W .* transpose(d.EOF) * cov(d.emb) * d.EOF))
end

"""
#save the results to jld2
"""
function save_results(d::matrixholder,name::String)
    jldsave("$name.jld2",
    EOF = Matrix(d.EOF),PC = Matrix(d.PC),RC = Matrix(d.RC),lambda = d.lambda,eps = d.eps, signal = d.signal)
end

"""
Lowpass filtering the signals : doing 7 for 7/a for now

savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/"*"fluxdata_raw.jld2"
wholedata = SharedArray{Float32}(load(savedirname)["data"])
wholedata_filtered = zeros(Float32,5114,18,16,3)
for i in 1:18
    for j in 1:16
        wholedata_filtered[:,i,j,1] = lowpass_filter(wholedata[:,i,j], 3)
        wholedata_filtered[:,i,j,2] = lowpass_filter(wholedata[:,i,j], 4)
        wholedata_filtered[:,i,j,3] = lowpass_filter(wholedata[:,i,j], 6)

    end
end

jldsave("/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_filtered.jld2",
       data_3a=wholedata_filtered[:,:,:,1],
       data_4a=wholedata_filtered[:,:,:,2],
       data_6a=wholedata_filtered[:,:,:,3])
"""

function lowpass_filter(signal::Vector{T}, cutoff_frequency) where T<:Real
    # Remove the mean of the signal
    signal_mean = mean(signal)
    signal .-= signal_mean

    # Compute the Fourier transform of the signal
    fourier = fft(signal) |> fftshift

    # Compute the frequencies corresponding to the Fourier coefficients
    Ts = 1 / 365.25
    t = 0:Ts:(N-1)*Ts
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

    # Find the indices of the Fourier coefficients corresponding to frequencies above the cutoff
    inds = findall(x -> abs(x) > cutoff_frequency, freqs)

    # Zero out the Fourier coefficients corresponding to frequencies above the cutoff
    fourier[inds] .= 0

    # Compute the inverse Fourier transform of the modified Fourier transform
    backtransform = real.(ifft(fourier |> fftshift)) .+ signal_mean

    return real(backtransform)
end
