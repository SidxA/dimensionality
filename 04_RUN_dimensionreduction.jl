"""
export JULIA_NUM_THREADS=1
GKSwstype=nul /opt/julia-1.9.0/bin/julia --threads 64 iterate_blockdata.jl 
"""

include("/net/home/lschulz/dimensionality/01_basic_functionality.jl")


#big parameters
N = 5114
k = 48
"""
running raw at the moment
"""

#fluxnetfullset
savedirname = "/net/scratch/lschulz/data/time_series.jld2"
#choose data: unfiltered, 3a,4a,6a lowpass filters
wholedata = SharedArray{Float32}(load(savedirname)["data_raw"])
#woledata = SharedArray{Float32}(load(savedirname)["data_3a"])
#wholedata = SharedArray{Float32}(load(savedirname)["data_4a"])
#wholedata = SharedArray{Float32}(load(savedirname)["data_6a"])

# Retrieve the dimensions of the loaded data
N, spots, vars = size(wholedata)


"""
OUTDIR NEEDS TO HAVE /
"""
outdir = "/net/scratch/lschulz/data/dimensionreduction/raw/"

"""
put EVERYTHING variable TO KNOW in the package
"""
#loc_ind var_ind
struct package
    outdir::String  # Directory to store the output
    loc_ind::Int64  # Location index
    var_ind::Int64  # Variable index
    
    function package(outdir, loc_ind, var_ind)
        # Constructor function to create a new 'package' object
        return new(outdir, loc_ind, var_ind)
    end
end


function packaging(outdir)
    # Check if the output directory exists, create it if it doesn't
    if isdir(outdir) == false
        mkdir(outdir)
    end
    
    P = package[]  # Initialize an empty array of package objects
    
    # Iterate over each spot and variable to create package objects and add them to the array
    for spot in 1:spots, vari in 1:vars
        P = push!(P, package(outdir, spot, vari))
    end
    
    return P  # Return the array of package objects
end




#iterate the packages of all specified parameters

function doitall(d::matrixholder, wholedata::SharedArray{Float32}, p::package)
    outdir = p.outdir
    loc_ind = p.loc_ind
    var_ind = p.var_ind

    # Set the sample year manually
    sampleyear = 365

    # Process the data using raw method
    put_data_raw(d, wholedata[:, loc_ind, var_ind])

    # Compute epsilon and perform differential EOF analysis
    put_epsilon(d)
    put_EOF_diff(d)
    calculate(d)

    # Save the results of differential EOF analysis
    save_results(d, outdir * "diff_$(d.W)_$(loc_ind)_$(var_ind)_raw")

    # Perform singular spectrum analysis (SSA) on the data
    put_EOF_ssa(d)
    calculate(d)

    # Save the results of SSA
    save_results(d, outdir * "ssa_$(d.W)_$(loc_ind)_$(var_ind)_raw")
end


#need to manually put the sampleyear in here...
#Walist = [5,5+1/4,5+1/3,5+1/2,5+2/3,5+3/4,6,6+1/4,6+1/3,6+1/2,6+2/3,6+3/4,7]
#Walist = [15,16,17]

Threads.@threads for p in packaging(outdir)
    # Iterate over the package objects in parallel using threads
    for W in Int.(floor.(365.25 .* [1, 7]))
        # Iterate over the desired values of W
        d = matrixholder(N, W, k)  # Create a new matrixholder object with the specified parameters
        doitall(d, wholedata, p)  # Call the doitall function with the matrixholder, wholedata, and package objects
        println("done $(p.loc_ind) $(p.var_ind) $(W)")  # Print a message to indicate the completion of the iteration
    end
end
