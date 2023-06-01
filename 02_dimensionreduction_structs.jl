
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


function put_data_gSSA(d::matrixholder, signal::Vector{Float32}, sampleyear::Int64)
    # Apply centralizer and deseasonalize functions to the signal using gSSA method
    # and assign the result to the 'signal' field of the matrixholder object
    d.signal = centralizer(deseasonalize_gSSA(centralizer(signal), sampleyear::Int64))
    
    # Compute the lagged embedding of the signal using the 'centralized_embed_lag' function
    # and assign the result to the 'emb' field of the matrixholder object
    d.emb = centralized_embed_lag(d.signal, d.W)'
end




function put_data_lSSA(d::matrixholder, signal::Vector{Float32}, sampleyear::Int64)
    # Apply centralizer and deseasonalize functions to the signal using lSSA method
    # and assign the result to the 'signal' field of the matrixholder object
    d.signal = centralizer(deseasonalize_lSSA(centralizer(signal), sampleyear::Int64))
    
    # Compute the lagged embedding of the signal using the 'centralized_embed_lag' function
    # and assign the result to the 'emb' field of the matrixholder object
    d.emb = centralized_embed_lag(d.signal, d.W)'
end


function put_data_window(d::matrixholder, signal::Vector{Float32}, sampleyear::Int64)
    # Calculate the window length based on a fraction of the sample year
    windowlength = Int(floor(sampleyear / 6))
    
    # Calculate the overlap between consecutive windows based on a fraction of the sample year
    overlap = Int(floor(sampleyear / 12))
    
    # Apply centralizer and moving_window_deseason functions to the signal using the specified window length and overlap
    # and assign the result to the 'signal' field of the matrixholder object
    d.signal = centralizer(moving_window_deseason(centralizer(signal), windowlength, overlap))
    
    # Compute the lagged embedding of the signal using the 'centralized_embed_lag' function
    # and assign the result to the 'emb' field of the matrixholder object
    d.emb = centralized_embed_lag(d.signal, d.W)'
end


function put_data_raw(d::matrixholder, signal::Vector{Float32})
    # Apply centralizer to the signal and assign the result to the 'signal' field of the matrixholder object
    d.signal = centralizer(signal)
    
    # Compute the lagged embedding of the signal using the 'centralized_embed_lag' function
    # and assign the result to the 'emb' field of the matrixholder object
    d.emb = centralized_embed_lag(d.signal, d.W)'
end


function centralized_embed_lag(data::Vector{Float32}, W::Int64)
    # Define a nested function 'centralizer' to centralize the data
    function centralizer(data)
        m = mean(data)  # Calculate the mean of the data
        cv = std(data)  # Calculate the standard deviation of the data
        return (data .- m) ./ cv  # Subtract the mean and divide by the standard deviation
    end
    
    Y = []  # Initialize an empty array to store the centralized data
    
    # Iterate over the data to create lagged embedding
    for i = 1:length(data) - W + 1
        Y = append!(Y, centralizer(data[i:i+W-1]))  # Centralize the data within the current window and append to 'Y'
    end
    
    # Reshape the array 'Y' into a matrix with 'W' rows and (length(data) - W + 1) columns
    # and convert the elements to Float32 type
    return reshape(float.(Y), W, length(data) - W + 1)::Matrix{Float32}
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
    
    # Constructor function for the 'epsilon' struct
    function epsilon(W, stepnumber, sampling_size)
        # Define the minimum and maximum values for 'eps' calculation
        min_eps = Float32(10^0)
        max_eps = Float32(10^7)

        # Calculate 'eps' values on a logarithmic scale between 'min_eps' and 'max_eps'
        eps = 10 .^ range(log10(min_eps), log10(max_eps), length = stepnumber)
        
        # Initialize 'data_samples' as an uninitialized matrix of size (W, sampling_size)
        data_samples = Matrix{Float32}(undef, W, sampling_size)
        
        # Initialize 'L' as an uninitialized vector of length 'eps'
        L = Vector{Float32}(undef, length(eps))
        
        # Initialize 'Weight' as an uninitialized matrix of size (sampling_size, sampling_size)
        Weight = Matrix{Float32}(undef, sampling_size, sampling_size)

        return new(data_samples, eps, L, Weight)
    end
end


# iteration function that returns the median
function fit_epsilon(data::Matrix{Float32}, stepnumber, sampling_size, it_number, W)
    P = size(data)[2]  # Number of columns in the data matrix
    object = epsilon(W, stepnumber, sampling_size)  # Create an instance of the 'epsilon' struct
    fit_eps_L = Vector{Float32}(undef, it_number)  # Initialize a vector to store the fit epsilon values
    sample = rand(1:P, sampling_size, it_number)  # Generate random samples from the columns of the data matrix
    
    for t in 1:it_number
        object.data_samples = data[:, sample[:, t]]  # Assign the sampled data to 'data_samples' field of 'object'
        
        for (i, eps) in enumerate(object.eps)
            for i in 1:sampling_size, j in 1:sampling_size
                # Compute the weight between data samples using the norm and eps
                object.Weight[i, j] = exp(-norm(object.data_samples[:, i] - object.data_samples[:, j])^2 / eps)
            end
            object.L[i] = sum(object.Weight)  # Compute the sum of weights for each eps
        end
        
        p0 = ones(3)  # Initial parameter values for curve fitting
        model(eps, p) = p[1] .* atan.(eps .- p[2]) .+ p[3]  # Define the model function for curve fitting
        
        # Fit the model function to log10(eps) and log10(L) using curve_fit and obtain the parameters
        p_midpoint = coef(curve_fit(model, log10.(object.eps), log10.(object.L), p0))
        
        a = p_midpoint[2]  # Extract the parameter 'a' from the fitted parameters
        b = median(log10.(object.L[1:8]))  # Compute the median of the log10(L) values
        one_over_e = model(b + (a - b) / exp(1), p_midpoint)  # Compute the value of the model at log10(1/e)
        
        fit_eps_L[t] = 10^(one_over_e)  # Compute the fit epsilon value for the current iteration
    end
    
    return Int64(floor(median(fit_eps_L)))  # Return the median of the fit epsilon values as an Int64
end


#function that works on object
function put_epsilon(d::matrixholder)
    stepnumber = 32  # Number of steps for epsilon calculation
    sampling_size = 64  # Size of the samples used for fitting epsilon
    it_number = 4  # Number of iterations for fitting epsilon

    # Calculate and assign the epsilon value to the 'eps' field of the matrixholder object
    d.eps = fit_epsilon(d.emb, stepnumber, sampling_size, it_number, d.W)
end

"""
#functions that fit in the model
"""
#little struct to hold all the variable types for the fit?
struct diffholder
    f::DiffMap{Float32}
    
    # Constructor function for the 'diffholder' struct
    function diffholder(data::Matrix{Float32}, k, t, alpha, eps)
        # Create a new instance of the 'diffholder' struct with a DiffMap object 'f' as a field
        return new(fit(DiffMap, data, maxoutdim=k, t=1, α=alpha, ɛ=eps))
    end
end


#type-stable extraction
function gettheproj_diff(diff::diffholder)
    # Return the projection matrix and eigenvalues from the 'diffholder' object
    return Matrix(diff.f.proj')::Matrix{Float32}, diff.f.λ::Vector{Float32}
end


#function that works on the object                  # FIXED t=1 alpha = 1
function put_EOF_diff(d::matrixholder)
    t = 1  # Time delay parameter
    alpha = 1.0  # Alpha parameter
    eps = d.eps  # Epsilon value

    # Create a 'diffholder' object by fitting a DiffMap to the lagged embedding 'd.emb'
    diff = diffholder(d.emb, d.k, t, alpha, eps)
    
    # Retrieve the projection matrix and eigenvalues from the 'diffholder' object
    d.EOF, d.lambda = gettheproj_diff(diff)
end


struct ssaholder
    f::PCA{Float32}
    
    # Constructor function for the 'ssaholder' struct
    function ssaholder(data::Matrix{Float32}, k)
        # Create a new instance of the 'ssaholder' struct with a PCA object 'f' as a field
        return new(fit(PCA, data', maxoutdim=k, method=:auto, pratio=1))
    end
end


#type-stable extraction
function gettheproj_ssa(ssa::ssaholder)
    # Return the projection matrix and normalized principal variances from the 'ssaholder' object
    return Matrix(ssa.f.proj)::Matrix{Float32}, (ssa.f.prinvars ./ ssa.f.prinvars[1])::Vector{Float32}
end


#function that works on the object
function put_EOF_ssa(d::matrixholder)
    d.EOF = zeros(Float32, d.W, d.k)  # Initialize EOF matrix with zeros
    d.lambda = zeros(Float32, d.k)  # Initialize lambda vector with zeros

    # Create an 'ssaholder' object by fitting a PCA to the lagged embedding 'd.emb'
    ssa = ssaholder(d.emb, d.k)
    
    # Retrieve the projection matrix and normalized principal variances from the 'ssaholder' object
    EOF, lambda = gettheproj_ssa(ssa)

    # Assign the retrieved EOF matrix to the 'd.EOF' field of the 'matrixholder' object 'd'
    d.EOF[:, 1:size(EOF)[2]] = EOF
    
    # Assign the retrieved lambda vector to the 'd.lambda' field of the 'matrixholder' object 'd'
    d.lambda[1:size(lambda)[1]] = lambda
end


"""
#building PC,RC   
"""                  
function calculate(d::matrixholder)
    # Calculate the principal components (PC) by multiplying the lagged embedding 'd.emb' with the EOF matrix 'd.EOF'
    d.PC = d.emb * d.EOF
    
    # Calculate the reconstructed components (RC) by iterating over each PC and its corresponding EOF
    d.RC = hcat([reconstructor(d.PC[:, i], d.EOF[:, i], d.N, d.W) for i in 1:d.k]...)
    
    # Calculate the lambda values (normalized eigenvalues) by applying the formula
    d.lambda = Float32.(diag(1/d.W .* transpose(d.EOF) * cov(d.emb) * d.EOF))
end


"""
#save the results to jld2
"""
function save_results(d::matrixholder, name::String)
    # Save the results to a JLD2 file with the specified name
    jldsave("$name.jld2", EOF = Matrix(d.EOF), PC = Matrix(d.PC), RC = Matrix(d.RC), lambda = d.lambda, eps = d.eps, signal = d.signal)
end

