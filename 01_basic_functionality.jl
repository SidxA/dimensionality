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
