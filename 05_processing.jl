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

#old ones
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
    return Int64.(ssa_h), Int64.(nlsa_h), ssa_trends, nlsa_trends
end

#only counting complete pairs
function create_trend_tensors_v2(outdir)

        #individual time series analysis
    function local_parameters(spot,vari,outdir)

        normalizer(x) = x ./ maximum(x)

        function gauss(x, p)
            # Gaussian function with parameters x0, gamma, and sigma
            x0 = p[1]
            gamma = p[2]
            sigma = p[3]
            @. return gamma * exp.(-(x - x0)^2 / sigma^2 / 2)
        end
        
        function fit_gauss(yvalues, xvalues)
            # Fit the Gaussian function to the data points
            onebins = xvalues
            bins = yvalues
            p0 = ones(3) .* 5
            return coef(curve_fit(gauss, bins, onebins, p0))
        end
        
        function harmonic_gaussian_per_mode(mode_spec, freqstart_w, freqend_w, freqs_w)
            # Fit a Gaussian model to a specific mode's spectrum within the specified frequency range
            spec = mode_spec
            spec[1:freqstart_w] .= 0
            spec[freqend_w:end] .= 0
        
            try
                # Attempt to fit the Gaussian function to the spectrum
                return fit_gauss(freqs_w, spec)
            catch
                # Return [0, 0, 0] if fitting fails
                return [0, 0, 0]
            end
        end
        
        
        function harmonicity_gauss(gausslist,eof_spec,freqstart_w,freqs_w)
        
            # Initialize empty arrays to store the harmonic and mixed frequency components, as well as their respective frequencies and residual components
            li_harmonics = Int64[] # array to store the indices of harmonic frequency components
            li_mixed = Int64[] # array to store the indices of mixed frequency components
            li_h_freq = Float64[] # array to store the frequencies of harmonic frequency components
            li_m_freq = Float64[] # array to store the frequencies of mixed frequency components
            li_residual=Int64[] # array to store the indices of frequency components with high residual values
            
            # Loop through each frequency component in the spectral components matrix eof_spec
            for i in 1:k
                # Get the ith spectral component
                mode = eof_spec[:,i]
                
                # Set the values of the first freqstart_w elements of the mode array to zero
                mode[1:freqstart_w] .= 0
                
                # Get the parameters of the Gaussian for the ith frequency component
                freq, value,sigma = gausslist[i]
        
                # Compute the Gaussian peak using the Gaussian function gauss defined earlier and the freqs_w array
                peak = gauss(freqs_w,(freq,value,sigma))
                
                # Compute the residual by subtracting the Gaussian peak from the mode array
                residual = mode .- peak
        
                # Determine if the ith frequency component is harmonic, mixed, or has a high residual value
                if maximum(residual .+ 0.0)/threshold1 <= value &&  any(abs.((1:8) .- freq).<=threshold2)
                    li_harmonics = append!(li_harmonics,i)
                    li_h_freq = append!(li_h_freq,freq)
                elseif maximum(residual .+ 0.0)/threshold1 >= value && any(abs.((1:8) .- freq).<=threshold2)
                    li_mixed = append!(li_mixed,i)
                    li_m_freq = append!(li_m_freq,freq)
                elseif maximum(residual .+ 0.0)/threshold1 >= value
                    li_residual = append!(li_residual,i)
                else
                    #println("no peak")
                end
        
            end
        
            # Return the harmonic frequency component indices
            return li_harmonics#,li_mixed,li_h_freq,li_m_freq,li_residual
        end

        #depending on N
        # Time parameters
        Ts = 1 / 365.25  # Time step (in years)
        t0 = 0  # Initial time
        tmax = t0 + (N - 1) * Ts  # Maximum time
        t = t0:Ts:tmax  # Time vector

        # Frequency parameters
        freqs = fftfreq(length(t), 1.0 / Ts) |> fftshift  # Frequency values
        freqstart = findall(x -> x >= 1 / 12, freqs)[1]  # Index of the starting frequency
        freqend = findall(x -> x >= 6, freqs)[1]  # Index of the ending frequency
        freq_domain_N = freqs[freqstart:freqend]  # Frequency domain within the specified range

        # Year values
        years = ((1:N) ./ 365) .+ startyear  # Year values corresponding to each time step


        #depending on W
        # Time parameters for windowed data
        tw = t0:Ts:(t0 + (W - 1) * Ts)  # Time vector for the windowed data

        # Frequency parameters for windowed data
        freqs_w = fftfreq(length(tw), 1.0 / Ts) |> fftshift  # Frequency values for the windowed data
        freqstart_w = findall(x -> x >= 1 / 12, freqs_w)[1]  # Index of the starting frequency for the windowed data
        freqend_w = findall(x -> x >= 6, freqs_w)[1]  # Index of the ending frequency for the windowed data
        freq_domain_w = freqs_w[freqstart_w:freqend_w]  # Frequency domain within the specified range for the windowed data


        #read in ssa nlsa results
        # Create the filenames for loading the results
        Filename_ssa = outdir * join(["ssa", W, spot, vari, preproc], "_") * "jld2"  # Filename for SSA results
        Filename_nlsa = outdir * join(["diff", W, spot, vari, preproc], "_") * "jld2"  # Filename for NLSA results

        # Load the SSA and NLSA result files
        file_ssa = load(Filename_ssa)  # Load the SSA result file
        file_nlsa = load(Filename_nlsa)  # Load the NLSA result file


        signal = file_ssa["signal"]

        # Extract signal from the SSA result file
        signal = file_ssa["signal"]

        # Extract SSA results
        ssa_lambda = file_ssa["lambda"]
        ssa_indices = sortperm(ssa_lambda, rev=true)
        ssa_Eof = file_ssa["EOF"][:, ssa_indices]
        ssa_PC = file_ssa["PC"][:, ssa_indices]
        ssa_RC = file_ssa["RC"][:, ssa_indices]
        ssa_lambda = ssa_lambda[ssa_indices]
        ssa_cap_var = ssa_lambda
        ssa_rec = ssa_RC

        # Extract NLSA results
        nlsa_lambda = file_nlsa["lambda"]
        nlsa_indices = sortperm(nlsa_lambda, rev=true)
        nlsa_Eof = file_nlsa["EOF"][:, nlsa_indices]
        nlsa_PC = file_nlsa["PC"][:, nlsa_indices]
        nlsa_RC = file_nlsa["RC"][:, nlsa_indices]
        nlsa_lambda = nlsa_lambda[nlsa_indices]
        nlsa_cap_var = nlsa_lambda

        nlsa_rec = nlsa_RC
        nlsa_eps = file_nlsa["eps"]


        #spectrum of signal and RC
        spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )
        spec_ssa_rc = (abs.(fft(ssa_rec) |> fftshift)[freqstart:freqend] |> normalizer )
        spec_nlsa_rc = (abs.(fft(nlsa_rec) |> fftshift)[freqstart:freqend] |> normalizer )

        #spectrum of individual eof
        spec_ssa_eof = hcat([abs.(fft(ssa_Eof[:,i]) |> fftshift)|> normalizer for i in 1:k]...) #[freqstart_w:freqend_w] 
        spec_nlsa_eof = hcat([abs.(fft(nlsa_Eof[:,i]) |> fftshift)|> normalizer for i in 1:k]...) #[freqstart_w:freqend_w] 

        #gaussian tables
        gaussian_ssa = [harmonic_gaussian_per_mode(spec_ssa_eof[:,i],freqstart_w,freqend_w,freqs_w) for i in 1:k]
        gaussian_nlsa = [harmonic_gaussian_per_mode(spec_nlsa_eof[:,i],freqstart_w,freqend_w,freqs_w) for i in 1:k]

        #harmonic indices
        li_harmonics_ssa = harmonicity_gauss(gaussian_ssa,spec_ssa_eof,freqstart_w,freqs_w)
        li_harmonics_nlsa = harmonicity_gauss(gaussian_nlsa,spec_nlsa_eof,freqstart_w,freqs_w)

        #seasonality behavior
        ssa_trend_harm = sum(ssa_RC[:,li_harmonics_ssa],dims=2)[:]
        nlsa_trend_harm = sum(nlsa_RC[:,li_harmonics_nlsa],dims=2)[:]

        #captured frequencies
        freq_ssa = [round(gaussian_ssa[i][1],digits=1) for i in li_harmonics_ssa]
        freq_nlsa = [round(gaussian_nlsa[i][1],digits=1) for i in li_harmonics_nlsa]

        #captured variance
        ssa_harm_var = round.(ssa_lambda[li_harmonics_ssa],digits=3)
        nlsa_harm_var = round.(nlsa_lambda[li_harmonics_nlsa],digits=3)

        #spectra of the seasonality and the residuals
        spec_ssa = (abs.(fft(ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer)
        spec_res_ssa = (abs.(fft(signal .- ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer)
        spec_nlsa = (abs.(fft(nlsa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer)
        spec_res_nlsa = (abs.(fft(signal .- nlsa_trend_harm) |> fftshift)[freqstart:freqend] |> normalizer )

        #varname
        varname = variables_names[vari]
        igbpclass = IGBP_list[spot]
        return [
            spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa
        ]
    end

    function fully_resolved(freq_list)
        b_ones = length(findall(x->abs(x-1.0)<=0.15,freq_list)) >= 2
        b_twos = length(findall(x->abs(x-2.0)<=0.15,freq_list)) >= 2
    
        #completely resolved
        if b_ones && b_twos
            return 3
        elseif b_ones
            return 2
        elseif length(freq_list) > 0
            return 1
        else
            return 0
        end
    end

    #needs to eat li_harmonics_nlsa
    function fully_resolved_trend_ind(freq_list,li_harmonics)
        b_ones = length(findall(x->abs(x-1.0)<=0.15,freq_list)) >= 2
        n_ones = findall(x->abs(x-1.0)<=0.15,freq_list)
        b_twos = length(findall(x->abs(x-2.0)<=0.15,freq_list)) >= 2
        n_twos = findall(x->abs(x-2.0)<=0.15,freq_list)
        b_threes = length(findall(x->abs(x-3.0)<=0.15,freq_list)) >= 2
        n_threes = findall(x->abs(x-3.0)<=0.15,freq_list)
        b_fours = length(findall(x->abs(x-4.0)<=0.15,freq_list)) >= 2
        n_fours = findall(x->abs(x-4.0)<=0.15,freq_list)
        b_fives = length(findall(x->abs(x-5.0)<=0.15,freq_list)) >= 2
        n_fives = findall(x->abs(x-5.0)<=0.15,freq_list)

        trend_ind = li_harmonics[cat(n_ones,n_twos,n_threes,n_fours,n_fives,dims=1)]
        return trend_ind
    end

    l_ssa = zeros(9,7)
    l_nlsa = zeros(9,7)
    ssa_trends = zeros(9, 7, N)
    nlsa_trends = zeros(9, 7, N)
    ssa_trends_pure = zeros(9, 7, N)
    nlsa_trends_pure = zeros(9, 7, N)

    for (i,spot) = enumerate(spots),(j,vari) = enumerate(vars)
        p = local_parameters(spot,vari,outdir)
        spoti, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p
        l_ssa[i,j] = fully_resolved(freq_ssa)
        l_nlsa[i,j] = fully_resolved(freq_nlsa)
        ssa_trend_inds = fully_resolved_trend_ind(freq_ssa,li_harmonics_ssa)
        nlsa_trend_inds = fully_resolved_trend_ind(freq_nlsa,li_harmonics_nlsa)
        ssa_trends_pure[i,j,:] = sum(ssa_rec[:,ssa_trend_inds],dims=2)[:]
        nlsa_trends_pure[i,j,:] = sum(nlsa_rec[:,nlsa_trend_inds],dims=2)[:]
        ssa_trends[i,j,:] = ssa_trend_harm
        nlsa_trends[i,j,:] = nlsa_trend_harm
    end

    return Int64.(l_ssa), Int64.(l_nlsa), ssa_trends_pure, nlsa_trends_pure, ssa_trends, nlsa_trends
end

function calculate_f_trends() # for v2
    trends_raw = create_trend_tensors_v2(outdir_raw)
    trends_f4 = create_trend_tensors_v2(outdir_f4)
    trends_f6 = create_trend_tensors_v2(outdir_f6)
    jldsave("/net/scratch/lschulz/data/trends.jld2",
        trends_raw = trends_raw,
        trends_f4 = trends_f4,
        trends_f6 = trends_f6,
    )
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
