include("/net/home/lschulz/dimensionality/01_basic_functionality.jl")

#global parameters
W = 2556
N = 5114
startyear = 2007
k = 48
preproc = "raw."
threshold1   = 0.15 #the harmonicity estimation
threshold2 = 0.15 #the harmonicity estimation

data = load("/net/scratch/lschulz/data/time_series.jld2")  # Load data from the JLD2 file

data_raw = data["data_raw"]  # unprocessed time series
data_f3 = data["data_f3"]  # lowpass filter at 3/a
data_f4 = data["data_f4"]  # lowpass filter at 4/a
data_f6 = data["data_f6"]  # lowpass filter at 6/a

spotslist = data["spotslist"]  # names of locations
IGBP_list = data["IGBP_list"]  # types of ecosystems

variables_names = data["variables_names"]  # names of variables
variables_original_names = data["variables_original_names"]  # fluxnet variable names

flags = data["flags"]  # quality flags
flag_variables = data["flag_variables"]  # names of the variables the quality flags are for


#clumsy coordinate trafo for indices
spots = mask_IGBP(IGBP_list)[1]
vars,varnames = mask_vari(variables_names)
IGBP_reduced = IGBP_list[spots]

function create_trend_tensors(N, W)
    # Initialize output matrices with zeros
    ssa_h = zeros(n_spots, n_vars)
    nlsa_h = zeros(n_spots, n_vars)
    ssa_trends = zeros(n_spots, n_vars, N)
    nlsa_trends = zeros(n_spots, n_vars, N)

    for spot in 1:n_spots
        for vari in 1:n_vars
            try
                # Get the parameters for the spot and variable
                p = local_parameters(W, vari, spot, N, startyear)
                # p = rescale_local_parameters(p)

                # Extract the necessary parameters from the parameter tuple
                spot, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p

                # Get the length of the harmonic vectors
                ssa_h[spot, vari] = length(p[end-5])
                nlsa_h[spot, vari] = length(p[end-4])

                # Get the trend vectors
                ssa_trends[spot, vari, :] = ssa_trend_harm
                nlsa_trends[spot, vari, :] = nlsa_trend_harm
            catch e
                # If an error occurs, do nothing
                # print("error")
                Nothing
            end
        end
    end
    return ssa_h, nlsa_h, ssa_trends, nlsa_trends
end

function regularity(signal)
    Ts = 1 / 365.25
    t0 = 0
    tmax = t0 + (N-1) * Ts
    t = t0:Ts:tmax

    # Compute the frequencies of the Fourier transform
    freqs = fftfreq(length(t), 1.0 / Ts) |> fftshift

    # Find the indices of the frequency range of interest
    freqstart = findall(x -> x >= 1/12, freqs)[1]
    freqend = findall(x -> x >= 6, freqs)[1]

    # Extract the frequency domain of interest
    freq_domain_N = freqs[freqstart:freqend]

    # Generate the years corresponding to the time series
    years = ((1:N) ./ 365) .+ startyear

    # Compute the Fourier transform of the signal
    spec_signal = abs.(fft(signal) |> fftshift)

    # Normalize the spectral power
    spec_signal /= norm(spec_signal)

    # Extract specific frequency components of interest
    flist = [findall(x -> (x - i)^2 < 0.001, freqs)[1] for i in 1:4]

    # Return the specific frequency components and the sum of the remaining components
    return spec_signal[flist], sum(spec_signal[freqend:end])
end

function long_deviation(vec) #with a/2
    st = std(vec)
    me = mean(vec)
    window_length = 182

    # Embed the lagged vector in a matrix
    Ar = embed_lag(Float32.(vec), window_length)

    for i in 1:size(Ar)[2]
        win = Ar[:, i]
        if all(win .< me)
            return true
        end
    end

    return false
end

function intra_regularity(signal)
    Ts = 1 / 365.25
    t0 = 0
    tmax = t0 + (N-1) * Ts
    t = t0:Ts:tmax

    # Compute the frequencies of the Fourier transform
    freqs = fftfreq(length(t), 1.0 / Ts) |> fftshift

    freqstart = findall(x -> x >= 1/12, freqs)[1]
    freqend = findall(x -> x >= 6, freqs)[1]

    # Compute the Fourier transform of the signal
    spec_signal = abs.(fft(signal) |> fftshift)

    h_freqs = [findall(x -> (x - i)^2 < 0.001, freqs)[1] for i in 1:4]
    harmonic_power = sum(spec_signal[h_freqs])
    noise_power = sum(spec_signal[freqend:end])
    power_sum = sum(spec_signal)

    return harmonic_power / power_sum, noise_power / power_sum
end

# Signal tensor-based analytics
function tensor_based_analytics(tensor)
    N, n_spots, n_varis = size(tensor)

    tensors = [intra_regularity(tensor[:, i, j]) for i = 1:n_spots, j = 1:n_varis] 
    harmonic_power = [tensors[i, j][1] for i in 1:n_spots, j = 1:n_varis]
    noise_power = [tensors[i, j][2] for i in 1:n_spots, j = 1:n_varis]

    return harmonic_power, noise_power
end

# Projection binning with few bins: for different heatmap symbols
#does not compare different filtered datasets with each other
function projection_binning_intra(matrices)
    function doit(matrix)
        n_spots, n_varis = size(matrix)
        bins = 3
        maxi = maximum(matrix) * 1.01
        mini = minimum(matrix)
        binsize = (maxi - mini) / bins
        bin_thresholds = [mini + i * binsize for i in 0:bins]

        binning = zeros(Int64, n_spots, n_varis)
        for i = 1:n_spots, j = 1:n_varis
            w = fit(Histogram, [matrix[i, j]], bin_thresholds).weights
            binning[i, j] = sum(w .* (1:bins))
        end

        return Int64.(binning)
    end

    return doit.(matrices)
end

# Intra-comparison is not suitable for filtered noise
# Build comparison across all filters
# create relative bins based on all variables, all locations, at all filters
function projection_binning_inter(matrices)
    n_spots, n_varis = size(matrices[1])
    n_matrices = size(matrices)[1]
    bins = 3
    maxi = maximum(maximum.(matrices)) * 1.01
    mini = minimum(minimum.(matrices))
    binsize = (maxi - mini) / bins
    bin_thresholds = [mini + i * binsize for i in 0:bins]

    binning = zeros(Int64, n_matrices, n_spots, n_varis)
    for m = 1:n_matrices, i = 1:n_spots, j = 1:n_varis
        w = fit(Histogram, [matrices[m][i, j]], bin_thresholds).weights
        binning[m, i, j] = sum(w .* (1:bins))
    end

    return [binning[m, :, :] for m = 1:n_matrices]
end

#required for sample entropy
function calculate_match_matrix(L, m, r)
    N = length(L)
    # Create a boolean matrix to store match information
    match_matrix = falses(N - m + 1, N - m + 1)
    for i in 1:N - m
        for j in i + 1:N - m + 1
            # Check if the absolute difference between L[i] and L[j] is within the threshold r
            if abs(L[i] - L[j]) ≤ r
                # If it satisfies the condition, mark it as a match
                match_matrix[i, j] = true
            else
                # Break the loop as further elements won't satisfy the condition
                break
            end
        end
    end
    return match_matrix
end

#required for sample entropy
function rle_encode(match_matrix)
    N = size(match_matrix, 1)
    encoded = Vector{Int}(undef, N)
    for i in 1:N
        encoded[i] = sum(match_matrix[i, i+1:end])
    end
    return encoded
end


#sample entropy
function sampen3(L)
    m = 2
    # Set the tolerance r as twice the standard deviation of the input sequence L
    r = 0.2 * std(L)
    N = length(L)
    B = 0
    A = 0

    # Calculate distance matrix for B
    match_matrix_B = calculate_match_matrix(L, m, r)
    encoded_B = rle_encode(match_matrix_B)
    # Sum the encoded values to get the count of B matches
    B = sum(encoded_B)

    # Calculate distance matrix for A
    m += 1
    match_matrix_A = calculate_match_matrix(L, m, r)
    encoded_A = rle_encode(match_matrix_A)
    # Sum the encoded values to get the count of A matches
    A = sum(encoded_A)

    # Handle division by zero to avoid errors
    if A == 0 || B == 0
        return NaN
    end

    # Calculate and return the Sample Entropy (SampEn)
    return -log(A / B)
end



function create_data_characteristics()
    # Compute entropy for each data set
    raw_entropy = [sampen3(data_raw[:,i,j]) for i = 1:18, j = 1:16]
    f4_entropy = [sampen3(data_f4[:,i,j]) for i = 1:18, j = 1:16]
    f6_entropy = [sampen3(data_f6[:,i,j]) for i = 1:18, j = 1:16]


    # Compute tensor-based analytics for raw data
    raw_harm_p, raw_noise_p = tensor_based_analytics(data_raw)

    # Compute tensor-based analytics for f4 data
    f4_harm_p, f4_noise_p = tensor_based_analytics(data_f4)

    # Compute tensor-based analytics for f6 data
    f6_harm_p, f6_noise_p = tensor_based_analytics(data_f6)

    # Extract noise, harmonic, and entropy values for each data set
    noises = [raw_noise_p[spots, vars], f4_noise_p[spots, vars], f6_noise_p[spots, vars]]
    harms = [raw_harm_p[spots, vars], f4_harm_p[spots, vars], f6_harm_p[spots, vars]]
    entropies = [raw_entropy[spots, vars], f4_entropy[spots, vars], f6_entropy[spots, vars]]

    # Scale the noise, harmonic, and entropy values using projection binning
    noises_scaled = projection_binning_intra(noises)
    harms_scaled = projection_binning_inter(harms)
    entropies_scaled = projection_binning_inter(entropies)

    artifacts = [long_deviation(flags[:,i,j]) for i = spots, j = [6,7,5,2,1,9,8]]
    # Save the computed data characteristics to a file
    jldsave("/net/scratch/lschulz/data/data_characteristics.jld2",
        f4_harm_p = f4_harm_p,
        f6_harm_p = f6_harm_p,
        raw_harm_p = raw_harm_p,
        f4_noise_p = f4_noise_p,
        f6_noise_p = f6_noise_p,
        raw_noise_p = raw_noise_p,
        raw_entropy = raw_entropy,
        f4_entropy = f4_entropy,
        f6_entropy = f6_entropy,
        raw_harm_p_scaled = harms_scaled[1],
        f4_harm_p_scaled = harms_scaled[2],
        f6_harm_p_scaled = harms_scaled[3],
        raw_noise_p_scaled = noises_scaled[1],
        f4_noise_p_scaled = noises_scaled[2],
        f6_noise_p_scaled = noises_scaled[3],
        raw_entropy_scaled = entropies_scaled[1],
        f4_entropy_scaled = entropies_scaled[2],
        f6_entropy_scaled = entropies_scaled[3],
        artifacts = artifacts,
    )
end
