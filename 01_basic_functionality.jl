using Pkg

#some weird error with datetime format, no time to set up new environment
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


#indices and names of main variables: ["GPP","RECO","NEE","SW_IN","TS","SWC"]
function mask_vari(variables_names)
    x = Int64[]  # Initialize an empty array to store the indices of the desired variables

    # Find the indices of specific variable names and append them to the array x
    x = append!(x, findall(x -> x == "GPP_DAY_1", variables_names))
    x = append!(x, findall(x -> x == "RECO_NIGHT_1", variables_names))
    x = append!(x, findall(x -> x == "NEE", variables_names))
    x = append!(x, findall(x -> x == "SW_IN", variables_names))
    x = append!(x, findall(x -> x == "TA", variables_names))
    x = append!(x, findall(x -> x == "TS", variables_names))
    x = append!(x, findall(x -> x == "SWC", variables_names))

    return x, ["GPP", "RECO", "NEE", "SW_IN", "TA", "TS", "SWC"]  # Return the array of indices and a new list of abbreviated variable names
end


#indices of main spots that are of forest ecosystem ([1]) and grass ecosystem ([2])
function mask_IGBP(IGBP_list)
    enf = findall(x -> x == "ENF", IGBP_list)  # Find the indices of "ENF" category in the IGBP_list
    mf = findall(x -> x == "MF", IGBP_list)  # Find the indices of "MF" 
    dbf = findall(x -> x == "DBF", IGBP_list)  # Find the indices of "DBF" 
    shr = findall(x -> x == "SHR", IGBP_list)  # Find the indices of "SHR" 
    cro = findall(x -> x == "CRO", IGBP_list)  # Find the indices of "CRO" 
    gra = findall(x -> x == "GRA", IGBP_list)  # Find the indices of "GRA" 
    osh = findall(x -> x == "OSH", IGBP_list)  # Find the indices of "OSH" 

    forest = append!(enf, mf, dbf)  # Combine the indices of forest categories into the forest array
    grass = append!(shr, gra, osh, cro)  # Combine the indices of grassland categories into the grass array

    return forest, grass  # Return the arrays of indices corresponding to forest and grass categories
end


#individual time series analysis based on custom saving directory
function extract_from_directory(saving_directory_ssa,saving_directory_nlsa)

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
    Filename_ssa = saving_directory_ssa * ".jld2"  # Filename for SSA results
    Filename_nlsa = saving_directory_nlsa * ".jld2"  # Filename for NLSA results

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
    vari = 1
    varname = "GPP"
    igbpclass = "IGBP"
    spot = 1
    return [
        spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa
    ]
end

#only counting complete pairs
function harmonic_structure(p_filtered,p_unfiltered)

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


    ssa_trends = zeros(2, N)
    nlsa_trends = zeros(2, N)
    ssa_harmonics = zeros(2, N)
    nlsa_harmonics = zeros(2, N)

    for (i,p) in enumerate([p_filtered, p_unfiltered])
        spoti, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p

        ssa_trend_inds = fully_resolved_trend_ind(freq_ssa,li_harmonics_ssa)
        nlsa_trend_inds = fully_resolved_trend_ind(freq_nlsa,li_harmonics_nlsa)
        ssa_trends_pure[i,:] = sum(ssa_rec[:,ssa_trend_inds],dims=2)[:]
        nlsa_trends_pure[i,:] = sum(nlsa_rec[:,nlsa_trend_inds],dims=2)[:]
        ssa_trends[i,:] = ssa_trend_harm
        nlsa_trends[i,:] = nlsa_trend_harm
    end

    return Int64.(l_ssa), Int64.(l_nlsa), ssa_trends_pure, nlsa_trends_pure, ssa_trends, nlsa_trends
end