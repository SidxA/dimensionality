#required functionality
include("/net/home/lschulz/dimensionality/fundamentals.jl")
using CairoMakie

#global parameters
W = 2556
N = 5114
startyear = 2007
k = 48
preproc = "raw."
threshold1   = 0.15 #the harmonicity estimation
threshold2 = 0.15 #the harmonicity estimation

#plotting settings
color_ssa = "darkgreen"
color_nlsa = "purple"
color_signal = "grey50"
lw = 3
lw_s = 1
ms = 5
fontsize = 22

#data
data = load("/net/scratch/lschulz/data/time_series.jld2")
data_raw = data["data_raw"]
data_f3 = data["data_f3" ]
data_f4 = data["data_f4"]
data_f6 = data["data_f6"]
spotslist = data["spotslist" ]
IGBP_list = data["IGBP_list"]
variables_names = data["variables_names"]
variables_original_names = data["variables_original_names"]
flags = data["flags"]
flag_variables = data["flag_variables"]

#individual dimension-reduction results
outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass4/"

"""
the chosen outdir needs to match the data_tensor choice in the functions!
"""

#seasonal cycle
cycles = load("/net/scratch/lschulz/data/seasonal_cycle.jld2")
ssa_h_raw = cycles["ssa_h_raw"]
ssa_h_3 = cycles["ssa_h_3"]
ssa_h_4 = cycles["ssa_h_4"]
ssa_h_6 = cycles["ssa_h_6"]
nlsa_h_raw = cycles["nlsa_h_raw"]
nlsa_h_3 = cycles["nlsa_h_3"]
nlsa_h_4 = cycles["nlsa_h_4"]
nlsa_h_6 = cycles["nlsa_h_6"]
ssa_trends_raw = cycles["ssa_trends_raw"]
ssa_trends_3 = cycles["ssa_trends_3"]
ssa_trends_4 = cycles["ssa_trends_4"]
ssa_trends_6 = cycles["ssa_trends_6"]
nlsa_trends_raw = cycles["nlsa_trends_raw"]
nlsa_trends_3 = cycles["nlsa_trends_3"]
nlsa_trends_4 = cycles["nlsa_trends_4"]
nlsa_trends_6 = cycles["nlsa_trends_6"]


#indices and names of main variables: ["GPP","RECO","NEE","SW_IN","TS","SWC"]
function mask_vari(variables_names)
    x = Int64[]
    x = append!(x,findall(x->x=="GPP_DAY_1",variables_names))
    x = append!(x,findall(x->x=="RECO_NIGHT_1",variables_names))
    x = append!(x,findall(x->x=="NEE",variables_names))
    x = append!(x,findall(x->x=="SW_IN",variables_names))
    x = append!(x,findall(x->x=="TS",variables_names))
    x = append!(x,findall(x->x=="SWC",variables_names))
    return x,["GPP","RECO","NEE","SW_IN","TS","SWC"]
end

#indices of main spots that are of forest ecosystem ([1]) and grass ecosystem ([2])
function mask_IGBP(IGBP_list)
    enf = findall(x->x=="ENF",IGBP_list)
    mf = findall(x->x=="MF",IGBP_list)
    dbf = findall(x->x=="DBF",IGBP_list)
    shr = findall(x->x=="SHR",IGBP_list)
    cro = findall(x->x=="CRO",IGBP_list)
    gra = findall(x->x=="GRA",IGBP_list)
    osh = findall(x->x=="OSH",IGBP_list)
    forest = append!(enf,mf,dbf)
    grass = append!(shr,gra,osh,cro)
    return forest,grass
end

#clumsy coordinate trafo for indices
spots = mask_IGBP(IGBP_list)[1]
vars,varnames = mask_vari(variables_names)
IGBP_reduced = IGBP_list[spots]

#data characteristics tensor
characteristics = load("/net/scratch/lschulz/data/data_characteristics.jld2")
f3_harm_p = characteristics["f3_harm_p"]
f4_harm_p = characteristics["f4_harm_p"]
f6_harm_p = characteristics["f6_harm_p"]
raw_harm_p = characteristics["raw_harm_p"]
f3_noise_p = characteristics["f3_noise_p"]
f4_noise_p = characteristics["f4_noise_p"]
f6_noise_p = characteristics["f6_noise_p"]
raw_noise_p = characteristics["raw_noise_p"]
raw_entropy = characteristics["raw_entropy"]
f3_entropy = characteristics["f3_entropy"]
f4_entropy = characteristics["f4_entropy"]
f6_entropy = characteristics["f6_entropy"]
raw_harm_p_scaled = characteristics["raw_harm_p_scaled"]
f3_harm_p_scaled = characteristics["f3_harm_p_scaled"]
f4_harm_p_scaled = characteristics["f4_harm_p_scaled"]
f6_harm_p_scaled = characteristics["f6_harm_p_scaled"]
raw_noise_p_scaled = characteristics["raw_noise_p_scaled"]
f3_noise_p_scaled = characteristics["f3_noise_p_scaled"]
f4_noise_p_scaled = characteristics["f4_noise_p_scaled"]
f6_noise_p_scaled = characteristics["f6_noise_p_scaled"]
raw_entropy_scaled = characteristics["raw_entropy_scaled"]
f3_entropy_scaled = characteristics["f3_entropy_scaled"]
f4_entropy_scaled = characteristics["f4_entropy_scaled"]
f6_entropy_scaled = characteristics["f6_entropy_scaled"]
artifacts = characteristics["artifacts"]

#individual time series analysis
function local_parameters(spot,vari,outdir)

    normalizer(x) = x./maximum(x)
    
    function gauss(x,p) # cauchy peak with two parameters
        x0 = p[1]
        gamma = p[2]
        sigma = p[3]
        #1/(sigma * sqrt(2*pi))
        @. return gamma * exp.(-(x-x0)^2/sigma^2/2)
    end
    
    function fit_gauss(yvalues,xvalues) # fit the cauchy peak to a vec of numbers
        onebins = xvalues
        bins = yvalues
        p0 = ones(3) .* 5
        return coef(curve_fit(gauss,bins,onebins,p0))
    end
    
    function harmonic_gaussian_per_mode(mode_spec,freqstart_w,freqend_w,freqs_w)
    
        spec = mode_spec
        spec[1:freqstart_w] .= 0
        spec[freqend_w:end] .= 0
    
        try
        return fit_gauss(freqs_w,spec)
        catch
            return [0,0,0]
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
    Ts = 1 / 365.25
    t0 = 0
    tmax = t0 + (N-1) * Ts
    t = t0:Ts:tmax
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    freq_domain_N = freqs[freqstart:freqend]
    years = ((1:N) ./ 365) .+ startyear

    #depending on W
    tw = t0:Ts:(t0+(W-1)*Ts)
    freqs_w = fftfreq(length(tw), 1.0/Ts) |> fftshift
    freqstart_w = findall(x->x>=1/12,freqs_w)[1]
    freqend_w = findall(x->x>=6,freqs_w)[1]
    freq_domain_w = freqs_w[freqstart_w:freqend_w]

    #read in ssa nlsa results
    Filename_ssa = outdir*join(["ssa", W, spot, vari, preproc], "_")*"jld2"
    Filename_nlsa = outdir*join(["diff", W, spot, vari, preproc], "_")*"jld2"
    file_ssa = load(Filename_ssa)
    file_nlsa = load(Filename_nlsa)

    signal = file_ssa["signal"]

    ssa_lambda = file_ssa["lambda"]
    ssa_indices = sortperm(ssa_lambda,rev=true)
    ssa_Eof = file_ssa["EOF"][:,ssa_indices]
    ssa_PC = file_ssa["PC"][:,ssa_indices]
    ssa_RC = file_ssa["RC"][:,ssa_indices]
    ssa_lambda = ssa_lambda[ssa_indices]
    #ssa_cap_var = sum(ssa_lambda)
    ssa_cap_var = ssa_lambda
    ssa_rec = sum(ssa_RC,dims=2)[:]

    nlsa_lambda = file_nlsa["lambda"]
    nlsa_indices = sortperm(nlsa_lambda,rev=true)
    nlsa_Eof = file_nlsa["EOF"][:,nlsa_indices]
    nlsa_PC = file_nlsa["PC"][:,nlsa_indices]
    nlsa_RC = file_nlsa["RC"][:,nlsa_indices]
    nlsa_lambda = nlsa_lambda[nlsa_indices]
    #nlsa_cap_var = sum(nlsa_lambda)
    nlsa_cap_var = nlsa_lambda

    nlsa_rec = sum(nlsa_RC,dims=2)[:]
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

#rescale time series to original scale in p, needs correct data-tensor, e.g. data_raw or data_f3
function rescale_local_parameters(p,data_tensor)

    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
    
    back_trafo(data,mean,std) = (data.*std).+mean

    m = mean(data_tensor[:,spot,vari])
    s = std(data_tensor[:,spot,vari])

    signal = back_trafo(signal,m,s)
    ssa_rec = back_trafo(ssa_rec,m,s)
    nlsa_rec = back_trafo(nlsa_rec,m,s)
    ssa_trend_harm = back_trafo(ssa_trend_harm,m,s)
    nlsa_trend_harm = back_trafo(nlsa_trend_harm,m,s)
    return [
        spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa
    ]
end
