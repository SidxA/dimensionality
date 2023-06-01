
include("/net/home/lschulz/dimensionality/fundamentals.jl")
using CairoMakie
dir = init_logging()

#global parameters
preproc = "raw." 
k= kappa = 48
yearsamples=365.25

#the harmonicity estimation
threshold1   = 0.15
threshold2 = 0.15

color_ssa = "darkgreen"
color_nlsa = "purple"
color_signal = "grey50"
lw = 3
ms = 5


odi = 4.5

#sites info

variables = ["GPP_DT_VUT_USTAR50","GPP_DT_VUT_MEAN","GPP_NT_VUT_USTAR50","GPP_NT_VUT_MEAN",
"RECO_DT_VUT_USTAR50","RECO_DT_VUT_MEAN","RECO_NT_VUT_USTAR50","RECO_NT_VUT_MEAN",
"NEE_VUT_MEAN","NEE_VUT_USTAR50_NIGHT","NEE_VUT_USTAR50_DAY",
"TA_F_DAY","SW_IN_F","LW_IN_F","SWC_F_MDS_1","TS_F_MDS_1"]

variables_names = ["GPP_DAY_1","GPP_DAY_2","GPP_NIGHT_1","GPP_NIGHT_2",
"RECO_DAY_1","RECO_DAY_2","RECO_NIGHT_1","RECO_NIGHT_2",
"NEE","NEE_NIGHT","NEE_DAY",
"TA","SW_IN","LW_IN","SWC","TS"]


#medium length
spotslist = [
    "BE-Lon",
    "BE-Vie",
    "CH-Dav",
    "CH-Fru",
    "CH-Lae",
    "CH-Oe2",
    "DE-Geb",
    "DE-Gri",
    "DE-Hai",
    "DE-Tha",
    "DK-Sor",
    "ES-LJu",
    "FI-Hyy",
    "FR-Aur",
    "FR-Lam",
    "IT-BCi",
    "IT-Lav",
    "RU-Fyo",
]

IGBP_list = [
    "CRO",
    "MF",
    "ENF",
    "GRA",
    "MF",
    "CRO",
    "CRO",
    "GRA",
    "DBF",
    "ENF",
    "DBF",
    "OSH",
    "ENF",
    "CRO",
    "CRO",
    "CRO",
    "ENF",
    "ENF",
]


year_ind =
    [
    1:365,
    366:730,
    731:1095,
    1096:1460,
    1462:1826, #89 1 offset
    1827:2191, 
    2192:2556,
    2557:2921,
    2923:3287, # 93 1 offset
    3288:3652,
    3653:4017,
    4018:4382,
    4384:4748, #97 1 offset
    4749:5113,
    #5114:5478,
    #5479:5843
]

lattitude = 

    [
    50.5516,
    50.3049,
    46.8153,
    47.1158,
    47.4783,
    47.2864,
    51.0997,
    50.9500,
    51.0792,
    50.9626,
    55.4859,
    36.9266,
    61.8474,
    43.549653,
    43.496437,
    40.5237,
    45.9562,
    56.4615,
]

longitude = 
    [
    4.7462,
    5.9981,
    9.8559,
    8.5378,
    8.3644,
    7.7337,
    10.9146,
    13.5126,
    10.4522,
    13.5651,
    11.6446,
    -2.7521,
    24.2948,
    1.106096,
    1.237878,
    14.9574,
    11.2813,
    32.9221,
]

elevation = 
    [
    167,
    493,
    1639,
    982,
    689,
    452,
    161.5,
    385,
    430,
    385,
    40,
    1600,
    181,
    250,
    181 ,
    20,
    1353,
    265,
]

outdir="/net/scratch/lschulz/fluxfullset_midwithnee_lowpass7/"
W = 2556
N = 5114
startyear = 2007


savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/"*"fluxdata_raw.jld2"
wholedata = SharedArray{Float32}(load(savedirname)["data"])

N,n_spots,n_vars = size(wholedata)

means = [mean(wholedata[:,spot,var]) for spot in 1:n_spots, var in 1:n_vars]
stds = [std(wholedata[:,spot,var]) for spot in 1:n_spots, var in 1:n_vars]


back_trafo(data,mean,std) = (data.*std).+mean

"""
functions for analysis
"""

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

"""
local parameters
"""



function local_parameters(W,vari,spot,N,startyear)

    normalizer(x) = x./maximum(x)
    
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


#rescale the signal, rc,  trend to the original scale

function rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
    signal = back_trafo(signal,means[spot,vari],stds[spot,vari])
    ssa_rec = back_trafo(ssa_rec,means[spot,vari],stds[spot,vari])
    nlsa_rec = back_trafo(nlsa_rec,means[spot,vari],stds[spot,vari])
    ssa_trend_harm = back_trafo(ssa_trend_harm,means[spot,vari],stds[spot,vari])
    nlsa_trend_harm = back_trafo(nlsa_trend_harm,means[spot,vari],stds[spot,vari])
    return [
        spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa
    ]
end

"""
plotting functionality
"""


# first one: create reconstruction picture
function plot_reconstruction(F,p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    #figure

    ax_time = Axis(F[1,1],xticks=Int.(floor.(years[1]:3:years[end])),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time (a)",
    ylabel=varname)

    ax_spec = Axis(F[2,1],yscale=log10,
    xlabel="frequency (/a)",
    ylabel="relative power",)

    scatter!(ax_time,years,signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_time,years,signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_time,years,ssa_rec,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_time,years,nlsa_rec,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain_N,spec_signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_spec,freq_domain_N,spec_ssa_rc,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain_N,spec_nlsa_rc,color=color_nlsa,linewidth=lw,label="NLSA")

    axislegend(ax_spec)

    text!(
        ax_time, 0, 1,
        text = "signal and reconstructions", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    text!(
        ax_spec, 0, 1,
        text = "spectra of signal and reconstructions", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    text!(
        ax_time, 0.95, 0.95,
        text = "$(spotslist[spot])\n$(igbpclass)\n$(varname)", 
        align = (:right, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = 20
    )

    return F
end

#second one: create modes picture
function plot_modeshapes(F,p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    ax_modes = Axis(F[1,1:3],yticksvisible = false,
    yticklabelsvisible = false,xticks=1:Int(floor(W/365.25)),
    xlabel="time (a)",
    ylabel="individual modes")


    ax_freq = Axis(F[1,4:6],limits=(1/10, 7,-0.7,37),#,yscale=log10
    xticks=1:7,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency (/a)",
    ylabel="relative power")

    ax_labels = Axis(F[1,7])
    hidedecorations!(ax_labels)

    #modes 

    color_fit = "grey"
    for i = 1:16

        mode = ssa_Eof[:,i]
        Four = spec_ssa_eof[:,i]
        #Four ./= maximum(abs.(Four)) 

        years_modes = (1:length(mode)) ./ 365.25


        lines!(ax_modes,years_modes,mode .+ (i*0.07),
        color=color_ssa)
        lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_ssa[i]) .+ i,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i, #[freqstart:end]
        color=color_ssa)

        mode = nlsa_Eof[:,i] 
        Four = spec_nlsa_eof[:,i]


        lines!(ax_modes,years_modes,mode .+ ((18+i)*0.07),
        color=color_nlsa)
        lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_nlsa[i]) .+ i .+ 18,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i .+ 18, #[freqstart:end]
        color=color_nlsa)

    end

    #height and freq
    boxstring = "peaks \n\n f \t height \n"
    for k = 16:-1:1
        boxstring *= string(round(gaussian_nlsa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_nlsa[k][2],digits=1))
        boxstring *= "\n"
    end
    boxstring*= "\n\n"
    for k = 16:-1:1
        boxstring *= string(round(gaussian_ssa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_ssa[k][2],digits=1))
        boxstring *= "\n"
    end

    text!(
        ax_labels, 0, 1,
        text = boxstring, 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    text!(
        ax_modes, 0, 1,
        text = "mode shapes", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    text!(
        ax_freq, 0, 1,
        text = "mode spectra", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    return F
end

#third one: create seasonality picture
function plot_seasonality(F,p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    seasonalitystring = ""
    seasonalitystring *= "ssa f \t"*string(round.(freq_ssa,digits=1))*"\n"
    seasonalitystring *= "nlsa f \t"*string(round.(freq_nlsa,digits=1))*"\n"
    seasonalitystring *= "ssa var \t"*string(round.(ssa_harm_var,digits=1))*"\n"
    seasonalitystring *= "nlsa var \t"*string(round.(nlsa_harm_var,digits=1))*"\n"

    #the figure

    offset = 2.5
    lw = 3
    lw_s = 1
    ms = 4

    textax = Axis(F[1,1])
    hidedecorations!(textax)
    text!(
        textax, 0, 1,
        text = seasonalitystring, 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    ax_time = Axis(F[2:3,1],
    xticks = Int.(floor.(years[1]:3:years[end])),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel = "time (a)",
    ylabel = varname,
    )


    #signal background

    scatter!(ax_time,years,signal .+ offset,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms,label="signal")
    lines!(ax_time,years,signal .+ offset,linewidth=lw_s,
    color=color_signal)
    #scatter!(ax_time,years,signal .+ 12.0,linewidth=lw,
    #color=color_signal,marker=:x,markersize=ms,label="signal")

    #ssa
    lines!(ax_time,years,ssa_trend_harm.+ offset,linewidth=lw,
    color=color_ssa,linestyle=:solid,label="ssa")
    lines!(ax_time,years,signal .- ssa_trend_harm .+ 0,linewidth=lw_s,
    color=color_ssa,linestyle=:solid,label="ssa")

    #nlsa
    lines!(ax_time,years,nlsa_trend_harm.+ offset,linewidth=lw,
    color=color_nlsa,linestyle=:solid,label="nlsa")
    lines!(ax_time,years,signal .- nlsa_trend_harm,linewidth=lw_s,
    color=color_nlsa,linestyle=:solid,label="nlsa")

    #spectrum
    ax_spec = Axis(F[4:5,1],
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(7),
    yscale = log10,
    xlabel = "frequency (/a)",
    ylabel = "relative power")

    #signal
    #scatter!(ax_spec,freq_domain,spec_signal,linewidth=lw,
    #color=color_signal,marker=:x,markersize=ms)
    scatter!(ax_spec,freq_domain_N,spec_signal .*10^odi,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms,label="signal")
    #lines!(ax_spec,freq_domain,spec_signal,linewidth=lw_s,
    #color=color_signal,linestyle=:solid)
    lines!(ax_spec,freq_domain_N,spec_signal .* 10^odi,linewidth=lw_s,
    color=color_signal,linestyle=:solid)

    #ssa
    lines!(ax_spec,freq_domain_N,spec_ssa.*10^odi,linewidth=lw,
    color=color_ssa,linestyle=:solid)
    lines!(ax_spec,freq_domain_N,spec_res_ssa,linewidth=lw,
    color=color_ssa,linestyle=:solid)

    #nlsa
    lines!(ax_spec,freq_domain_N,spec_nlsa.*10^odi,linewidth=lw,
    color=color_nlsa,linestyle=:solid)
    lines!(ax_spec,freq_domain_N,spec_res_nlsa,linewidth=lw,
    color=color_nlsa,linestyle=:solid)



    text!(
        ax_time, 0, 1,
        text = "seasonal behavior", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )
    text!(
        ax_time, 0, 0,
        text = "residual behavior", 
        align = (:left, :bottom),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )
    text!(
        ax_spec, 0, 1,
        text = "seasonal spectrum", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )
    text!(
        ax_spec, 0, 0,
        text = "residual spectrum", 
        align = (:left, :bottom),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )
    return F
end

#fourth one: create annual comparison picture
function plot_annuality(F,p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p




    years_plot = Int.(floor.(years[1]:1:years[end-10]))

    color_signal_grad =  cgrad([:grey20,:grey80],length(years),categorical=true,rev=true)
    color_diff_grad = color_signal_grad
    color_ssa_grad =  cgrad([:darkgreen,:lightgreen],length(years),categorical=true,rev=true)
    color_nlsa_grad =  cgrad([:purple,:plum],length(years),categorical=true,rev=true)

    lw = 1


    #fund = hcat([trend_first[i] for i in year_ind]...)

    ssa_harm = hcat([ssa_trend_harm[i] for i in year_ind]...)
    nlsa_harm = hcat([nlsa_trend_harm[i] for i in year_ind]...)

    #ax_years = Axis(F[1:2, 1],
    #    yticks = ((1:Int(floor(N/365.25))) ./ 2 .-1,  string.(years_plot)),
    #    xticks = ([365/2, 400 + 365/2],
    #    ["SSA harm","NLSA harm"]),
    #    ylabel="years"
    #)
    #for i in 1:length(years_plot)
    #    #d = lines!(ax_years,1:365,fund[:,i].+(0.5*i),color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
    #    d = lines!(ax_years,(1:365),ssa_harm[:,i].+(0.5*i),color = color_ssa_grad[i],linewidth=lw)
    #    #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
    #    d = lines!(ax_years,400 .+ (1:365),nlsa_harm[:,i].+(0.5*i),color = color_nlsa_grad[i],linewidth=lw)
    #end

    ax_years_2 = Axis(F[1:2,1],
    xticks = (vcat((365/4)*(0:1:3),(365/4)*(0:1:3).+400,(365/4)*(0:1:3).+800),
        ["Jan","Apr","Jul","Oct","Jan","Apr","Jul","Oct","Jan","Apr","Jul","Oct"]),
        ylabel=varname
        )

    for i in 1:length(years_plot)
        #d = lines!(ax_years_2,1:365,fund[:,i],color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
        
        d = lines!(ax_years_2,(1:365),ssa_harm[:,i],color = color_ssa_grad[i],linewidth=lw)
        d = lines!(ax_years_2,400 .+ (1:365),nlsa_harm[:,i],color = color_nlsa_grad[i],linewidth=lw)

        d = lines!(ax_years_2,800 .+ (1:365),ssa_harm[:,i] .- nlsa_harm[:,i],color = color_diff_grad[i],linewidth=lw)

    end

    text!(
        ax_years_2, 0, 0.1,
        text = "SSA", 
        align = (:left, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    text!(
        ax_years_2, 0.5,0.1,
        text = "NLSA", 
        align = (:left, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    text!(
        ax_years_2, 1, 0.1,
        text = "SSA - NLSA", 
        align = (:right, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = 12
    )

    text!(
        ax_years_2, 0, 1,
        text = "change of seasonal behavior", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    ax_years_3_maxima = Axis(F[3,1],xticks=years_plot,ylabel=varname,xlabel="time (a)")

    yearcomparison = years_plot

    #fund_yearmax = [maximum(fund[:,i]) for i in 1:16]
    ssa_yearmax = [maximum(ssa_harm[:,i]) for i in 1:length(years_plot)]
    nlsa_yearmax = [maximum(nlsa_harm[:,i]) for i in 1:length(years_plot)]

    #lines!(ax_years_3_maxima,yearcomparison,fund_yearmax,color = color_signal,linewidth=lw)
    #scatter!(ax_years_3_maxima,yearcomparison,fund_yearmax,color = color_signal,linewidth=lw,markersize = ms)
    lines!(ax_years_3_maxima,yearcomparison,ssa_yearmax,color = color_ssa,linewidth=lw)
    scatter!(ax_years_3_maxima,yearcomparison,ssa_yearmax,color = color_ssa,linewidth=lw,markersize = ms)
    lines!(ax_years_3_maxima,yearcomparison,nlsa_yearmax,color = color_nlsa,linewidth=lw)
    scatter!(ax_years_3_maxima,yearcomparison,nlsa_yearmax,color = color_nlsa,linewidth=lw,markersize = ms)

    text!(
        ax_years_3_maxima, 0, 1,
        text = "change of maximum value", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    ax_years_4_max_positions = Axis(F[4,1],xticks=years_plot,ylabel="day of year",xlabel="time (a)")

    halfyear = Int(floor(365/2))
    hlines!(ax_years_4_max_positions,halfyear , xmax = 1, color = :grey)

    #fund_argmax = [argmax(fund[:,i]) for i in 1:16]
    ssa_argmax = [argmax(ssa_harm[:,i]) for i in 1:length(years_plot)]
    nlsa_argmax = [argmax(nlsa_harm[:,i]) for i in 1:length(years_plot)]

    #lines!(ax_years_4_max_positions,yearcomparison,fund_argmax,color = color_signal,linewidth=lw)
    #scatter!(ax_years_4_max_positions,yearcomparison,fund_argmax,color = color_signal,linewidth=lw,markersize = ms)
    lines!(ax_years_4_max_positions,yearcomparison,ssa_argmax,color = color_ssa,linewidth=lw)
    scatter!(ax_years_4_max_positions,yearcomparison,ssa_argmax,color = color_ssa,linewidth=lw,markersize = ms)
    lines!(ax_years_4_max_positions,yearcomparison,nlsa_argmax,color = color_nlsa,linewidth=lw)
    scatter!(ax_years_4_max_positions,yearcomparison,nlsa_argmax,color = color_nlsa,linewidth=lw,markersize = ms)

    text!(
        ax_years_4_max_positions, 0, 1,
        text = "change of maximum position", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 16
    )

    #Label(F[0, :], text = "changes in seasonal behavior for "*titlestring,
    #fontsize = 24)

    #trim!(f.layout)


    return F
end

#put them all together
function create_combined_plot(spot,savedirname,W,vari,N)

        p = local_parameters(W,vari,spot,N,startyear)
        p = rescale_local_parameters(p)
    
        F = Figure(resolution=(1600,1600))

        plot_reconstruction(F[1,1],p)
        plot_modeshapes(F[1,2],p)
        plot_seasonality(F[2,1],p)
        plot_annuality(F[2,2],p)
    
        !isdir(savedirname) ? mkpath(savedirname) : "go" 
    
        savedirstring = "spot$(spot)_W$(W)_$(variables_names[vari])"
        savedir = savedirname*savedirstring

        save(savedir*".png",F)
    
end

#iterate over the spots and variables
function plot_existing_data()
    for vari = 1:n_vars, spot = 1:n_spots
        try
            savedirname = dir*"runs/spot$spot/"
            create_combined_plot(spot,savedirname,W,vari,N)
            println("YES spot $spot variable $vari")
        catch
            println("NO spot $spot variable $vari")
        end
    end
end

"""
summarizing quantification + plotting heatmaps
"""

#how does the value of the yearly maximum change linearly
function maxima_linear_slope(p)
    # Extract input parameters
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    # Calculate maximum year from the length of the time series
    maxyear = Int(floor(N/365))

    # Create an array of years from 1 to the maximum year
    x = Float64[1:maxyear...]

    # Extract harmonic trends for ssa and nlsa modes for all selected years
    ssa_harm = hcat([ssa_trend_harm[i] for i in year_ind]...)
    nlsa_harm = hcat([nlsa_trend_harm[i] for i in year_ind]...)

    # Calculate the maximum value of each year for ssa and nlsa
    ssa_yearmax = [maximum(ssa_harm[:,i]) for i in 1:maxyear]
    nlsa_yearmax = [maximum(nlsa_harm[:,i]) for i in 1:maxyear]

    # Calculate the index of the maximum value of each year for ssa and nlsa
    ssa_argmax = [argmax(ssa_harm[:,i]) for i in 1:maxyear]
    nlsa_argmax = [argmax(nlsa_harm[:,i]) for i in 1:maxyear]

    # Perform linear regression on the year and ssa/nlsa maximum values
    ssa_yearmax_linreg = llsq(x,Float64.(ssa_yearmax))[1]
    nlsa_yearmax_linreg = llsq(x,Float64.(nlsa_yearmax))[1]

    # Perform linear regression on the year and ssa/nlsa maximum value indices
    ssa_argmax_linreg = llsq(x,Float64.(ssa_argmax))[1]
    nlsa_argmax_linreg = llsq(x,Float64.(nlsa_argmax))[1]

    # Return the linear regression slopes for ssa and nlsa year-max and argmax values
    return [ssa_yearmax_linreg,nlsa_yearmax_linreg,ssa_argmax_linreg,nlsa_argmax_linreg]
end

#tables ssa/nsla separate
#nr of harmonics
#lin slope of year max
#lin slope of year argmax


function create_valuetables(N,W)
    # Initialize output matrices with zeros
    ssa_harmonics = zeros(Int64,n_spots,n_vars)
    nlsa_harmonics = zeros(Int64,n_spots,n_vars)
    ssa_max_slope = zeros(n_spots,n_vars)
    nlsa_max_slope = zeros(n_spots,n_vars)
    ssa_argmax_slope = zeros(n_spots,n_vars)
    nlsa_argmax_slope = zeros(n_spots,n_vars)

    for spot in 1:n_spots
        for vari in 1:n_vars
            try
                # Get the parameters for the spot and variable
                p = local_parameters(W,vari,spot,N,startyear)
                
                # Get the length of the harmonic vectors
                ssa_h = length(p[end-5])
                nlsa_h = length(p[end-4])
                ssa_harmonics[spot,vari] = ssa_h
                nlsa_harmonics[spot,vari] = nlsa_h
                
                ssa_yearmax_linreg,nlsa_yearmax_linreg,ssa_argmax_linreg,nlsa_argmax_linreg = maxima_linear_slope(p)

                if ssa_h >= 2
                    ssa_max_slope[spot,vari] = ssa_yearmax_linreg
                    ssa_argmax_slope[spot,vari] = ssa_argmax_linreg
                end
                if nlsa_h >= 2
                    nlsa_max_slope[spot,vari] = nlsa_yearmax_linreg
                    nlsa_argmax_slope[spot,vari] = nlsa_argmax_linreg
                end

            catch e
                # If an error occurs, do nothing
                Nothing
            end
        end
    end


    # Return the output matrices
    return ssa_harmonics,nlsa_harmonics,ssa_max_slope,nlsa_max_slope,ssa_argmax_slope,nlsa_argmax_slope
end


function heatmap_table(f,valuetable,titlestring,c = :solar)
    stepnr = 100
    stepsize = (maximum(valuetable)-minimum(valuetable))/(stepnr) 
    bins = LinRange(minimum(valuetable)-stepsize,maximum(valuetable)+stepsize,stepnr+2)

    cgradient = cgrad(c, length(bins), categorical = true)
    belongstobin(x) = findall(bins .- stepsize .<= x .< bins .+ stepsize)[end]

    ax = Axis(f[1,1],
    aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*"\n".*IGBP_list),
    xticklabelrotation = pi/2,
    #xticklabelalign = (:left,:top),
    xticklabelsize = 12,
    yticks = (1:length(variables_names),variables_names),
    #yticklabelalign = (:right,:center),
    yticklabelsize = 12,
    subtitle = titlestring
    )



    hm = heatmap!(ax,valuetable,
    colormap = cgradient,
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    aspect_ratio = 1,
    grid = true,
    framestyle = :box,
    )

    Colorbar(f[1,2], hm)

    return f
end

"""
phase and amplitude investigations
"""

protophase(signal) =  atan.(imag(analyticsignal(Float64.(signal))),real(analyticsignal(Float64.(signal))))
amplitude(signal) = median([maximum(hcat([signal[i] for i in year_ind]...)[:,i]) for i in 1: Int(floor(N/365))])
argmax_position(signal) =  median([argmax(hcat([signal[i] for i in year_ind]...)[:,i]) for i in 1: Int(floor(N/365))])

function amplitude_diff(signal)
    max = [maximum(hcat([signal[i] for i in year_ind]...)[:,i]) for i in 1: Int(floor(N/365))]
    min = [minimum(hcat([signal[i] for i in year_ind]...)[:,i]) for i in 1: Int(floor(N/365))]
    return median(max.-min)
end

    

#create tensors for spot,vari : seasonal trends and n_harmonics
function create_trend_tensors(N,W)
    # Initialize output matrices with zeros
    ssa_h = zeros(n_spots,n_vars)
    nlsa_h = zeros(n_spots,n_vars)
    ssa_trends = zeros(n_spots,n_vars,N)
    nlsa_trends = zeros(n_spots,n_vars,N)

    for spot in 1:n_spots
        for vari in 1:n_vars
            try
                # Get the parameters for the spot and variable
                p = local_parameters(W,vari,spot,N,startyear)
                #p = rescale_local_parameters(p)
                spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

                # Get the length of the harmonic vectors
                ssa_h[spot,vari] = length(p[end-5])
                nlsa_h[spot,vari] = length(p[end-4])
                
                # Get the trend vectors
                ssa_trends[spot,vari,:] = ssa_trend_harm
                nlsa_trends[spot,vari,:] = nlsa_trend_harm

            catch e
                # If an error occurs, do nothing
                #print("error")
                Nothing
            end
        end
    end
    return ssa_h,nlsa_h,ssa_trends,nlsa_trends
end

#bins of protophase - relative sin
function phase_diff(phases)
    rel_sin = sin.(2 .*pi./365 .*(1:N))
    phasediff = phases .- rel_sin
    #phasediff[phasediff .> pi] .-= 2*pi
    #phasediff[phasediff .< -pi] .+= 2*pi
    phasediff[phasediff .< 0] .+= 2*pi
    stepsize = 0.01
    bins = (0:stepsize:2*pi)
    h =fit(Histogram,phasediff,bins)
    h = normalize(h, mode=:pdf)
    return bins[1:end-1],h.weights'
end

#fit cauchy to phase difference
function phase_offset(bins,weights)

    function cauchy(x,p) # fit cauchy / delta to a normalized phase difference
        x0 = p[1]
        gamma = p[2]
        @. return 1/(pi*gamma*(1+((x-x0)/gamma)^2))
    end

    p0 = ones(2)
    t = curve_fit(cauchy,bins[:],weights[:],p0)
    params = t.param
    return params[1],params[2],cauchy(bins,params)
end

function nlsa_overlay(M)
    M[nlsa_h .< 2] .= NaN
    return M
end

function all_phase_plots()
        
    ssa_h,nlsa_h,ssa_trends,nlsa_trends = create_trend_tensors(N,W)
    jldsave(dir*"trends.jld2",ssa_h = ssa_h,nlsa_h = nlsa_h,ssa_trends = ssa_trends,nlsa_trends = nlsa_trends)

    ssa_ht = [protophase(ssa_trends[i,j,:]) for i in 1:n_spots, j in 1:n_vars]
    nlsa_ht = [protophase(nlsa_trends[i,j,:]) for i in 1:n_spots, j in 1:n_vars]

    ssa_cauchyfit = [phase_offset(phase_diff(ssa_ht[i,j])...) for i in 1:n_spots, j in 1:n_vars]
    ssa_offset = [ssa_cauchyfit[i,j][1] for i in 1:n_spots, j in 1:n_vars]
    ssa_offset_strength = [ssa_cauchyfit[i,j][2] for i in 1:n_spots, j in 1:n_vars]
    ssa_offset_fit = [ssa_cauchyfit[i,j][3] for i in 1:n_spots, j in 1:n_vars]

    nlsa_cauchyfit = [phase_offset(phase_diff(nlsa_ht[i,j])...) for i in 1:n_spots, j in 1:n_vars]
    nlsa_offset = [nlsa_cauchyfit[i,j][1] for i in 1:n_spots, j in 1:n_vars]|> nlsa_overlay
    nlsa_offset_strength = [nlsa_cauchyfit[i,j][2] for i in 1:n_spots, j in 1:n_vars] |> nlsa_overlay
    nlsa_offset_fit = [nlsa_cauchyfit[i,j][3] for i in 1:n_spots, j in 1:n_vars]

    ssa_amplitude = [amplitude_diff(ssa_trends[i,j,:]) for i in 1:n_spots, j in 1:n_vars]
    nlsa_amplitude = [amplitude_diff(nlsa_trends[i,j,:]) for i in 1:n_spots, j in 1:n_vars] |> nlsa_overlay

    nlsa_h = nlsa_h |> nlsa_overlay

    ssa_argmax_pos = [argmax_position(ssa_trends[i,j,:]) for i in 1:n_spots, j in 1:n_vars]
    nlsa_argmax_pos = [argmax_position(nlsa_trends[i,j,:]) for i in 1:n_spots, j in 1:n_vars] |> nlsa_overlay

    #plotting

    save(dir*"ssa_harmonics.png",heatmap_table(Figure(),ssa_h,"ssa harmonics"))
    save(dir*"nlsa_harmonics.png",heatmap_table(Figure(),nlsa_h,"nlsa harmonics"))
    save(dir*"ssa_offset_strength.png",heatmap_table(Figure(),ssa_offset_strength,"ssa offset strength"))
    save(dir*"nlsa_offset_strength.png",heatmap_table(Figure(),nlsa_offset_strength,"nlsa offset strength"))
    save(dir*"ssa_amplitude.png",heatmap_table(Figure(),ssa_amplitude,"ssa amplitude"))
    save(dir*"nlsa_amplitude.png",heatmap_table(Figure(),nlsa_amplitude,"nlsa amplitude"))


    save(dir*"ssa_offset.png",heatmap_table(Figure(),ssa_offset,"ssa offset", :cyclic_protanopic_deuteranopic_bwyk_16_96_c31_n256))
    save(dir*"nlsa_offset.png",heatmap_table(Figure(),nlsa_offset,"nlsa offset", :cyclic_protanopic_deuteranopic_bwyk_16_96_c31_n256))

    save(dir*"ssa_argmax_pos.png",heatmap_table(Figure(),ssa_argmax_pos,"ssa argmax position",:cyclic_protanopic_deuteranopic_bwyk_16_96_c31_n256))
    save(dir*"nlsa_argmax_pos.png",heatmap_table(Figure(),nlsa_argmax_pos,"nlsa argmax position",:cyclic_protanopic_deuteranopic_bwyk_16_96_c31_n256))
end

"""
bring it all together
1) fundamental vs harmonic
"""

#function to deliver individual harmonics
# needs to be apllied to non-rescaled p
function individual_harmonics(p,rescale = false)
    #redundant double calculation of rc because it simplifies code
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
    l = embed_lag(signal,W)'
    ssa_trend_rc = hcat([reconstructor((l * ssa_Eof)[:,i],ssa_Eof[:,i],N,W) for i in li_harmonics_ssa]...)
    nlsa_trend_rc = hcat([reconstructor((l * nlsa_Eof)[:,i],nlsa_Eof[:,i],N,W) for i in li_harmonics_nlsa]...)
    if rescale == true
        m = means[spot,vari]
        s = stds[spot,vari]
        ssa_trend_rc = hcat([back_trafo(ssa_trend_rc[:,i],m,s) for i in 1:size(ssa_trend_rc,2)]...)
        nlsa_trend_rc = hcat([back_trafo(nlsa_trend_rc[:,i],m,s) for i in 1:size(nlsa_trend_rc,2)]...)
    end
    return ssa_trend_rc,nlsa_trend_rc
end


#needs to be put in not rescaled p
function plot_variable_dynamics_rescaled(F,p,fontsize,rel_offset)
    #2 plots of series similar to trend reconstruction
    ssa_trend_rc,nlsa_trend_rc = individual_harmonics(p,false)
    p = rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
    
    m = means[spot,vari]
    s = stds[spot,vari]

    ssa_fund = back_trafo(sum(ssa_trend_rc[:,1:2],dims=2)[:],m,s)
    nlsa_fund = back_trafo(sum(nlsa_trend_rc[:,1:2],dims=2)[:],m,s)

    seasonalitystring = ""
    seasonalitystring *= "ssa f \t"*string(round.(freq_ssa,digits=1))*"\n"
    seasonalitystring *= "nlsa f \t"*string(round.(freq_nlsa,digits=1))*"\n"
    seasonalitystring *= "ssa var \t"*string(round.(ssa_harm_var,digits=1))*"\n"
    seasonalitystring *= "nlsa var \t"*string(round.(nlsa_harm_var,digits=1))*"\n"

    #the figure

    offset = (maximum(signal)-minimum(signal)) * rel_offset
    lw = 3
    lw_s = 1
    ms = 4


    ax_time = Axis(F[1,1],
    xticks = Int.(floor.(years[1]:3:years[end])),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel = "time (a)",
    ylabel = varname,
    xlabelsize = fontsize,
    ylabelsize = fontsize,
    xticklabelsize = fontsize-4,
    yticklabelsize = fontsize-4,
    )


    """
    H
    """
    #signal background

    scatter!(ax_time,years,signal .+ offset,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms)
    li = lines!(ax_time,years,signal .+ offset,linewidth=lw_s,
    color=color_signal,label="signal")

    #ssa
    li = lines!(ax_time,years,ssa_trend_harm .+ offset,linewidth=lw,
    color=color_ssa,linestyle=:solid,label="SSA")

    #nlsa
    li = lines!(ax_time,years,nlsa_trend_harm .+ offset,linewidth=lw,
    color=color_nlsa,linestyle=:solid,label="NLSA")

    """
    F
    """
    scatter!(ax_time,years,signal,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms)
    lines!(ax_time,years,signal,linewidth=lw_s,
    color=color_signal)

    #ssa
    lines!(ax_time,years,ssa_fund,linewidth=lw,
    color=color_ssa,linestyle=:solid)


    #nlsa
    lines!(ax_time,years,nlsa_fund,linewidth=lw,
    color=color_nlsa,linestyle=:solid)



    text!(
        ax_time, 0, 1,
        text = "harmonics", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fontsize
    )
    text!(
        ax_time, 0, 0.02,
        text = "fundamental", 
        align = (:left, :bottom),
        offset = (4, -2),
        space = :relative,
        fontsize = fontsize
    )

    text!(
        ax_time, 0.95, 1,
        text = "$(spotslist[spot]) $(igbpclass) $(varname)", 
        align = (:right, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = fontsize
    )

    hideydecorations!(ax_time, label = false)
    hidespines!(ax_time, :t, :r)

    tightlimits!(ax_time,Left(),Right())

    axislegend(ax_time,
        tellheight = false,
        tellwidth = false,
        margin = (4, 4, 4, 4),
        halign = :right, valign = :bottom, orientation = :horizontal,
        labelsize = fontsize-2,
    )
    return F
end


"""
2) heatmap: SSA vs NLSA
#take the H from the saved files
"""

function heatmap_comparison(F,ssa_h,nlsa_h)

    nlsa_h[nlsa_h .< 2] .= NaN
    ssa_h[ssa_h .< 2] .= NaN

    fs = 10
    ax = Axis(F[1,1],
    aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    #xticklabelalign = (:left,:top),
    xticklabelsize = fs,
    yticks = (1:length(variables_names),variables_names),
    #yticklabelalign = (:right,:center),
    yticklabelsize = fs,
    subtitle = "ssa"
    )



    hm = heatmap!(ax,ssa_h,
    colormap = :solar,
    colorrange = (0,12),
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    aspect_ratio = 1,
    grid = true,
    framestyle = :box,

    )

    ax = Axis(F[1,2],
    aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    #xticklabelalign = (:left,:top),
    xticklabelsize = fs,
    yticks = (1:length(variables_names),variables_names),
    #yticklabelalign = (:right,:center),
    yticklabelsize = fs,
    subtitle = "nlsa"
    )

    hm = heatmap!(ax,nlsa_h,
    colormap = :solar,
    colorrange = (0,12),
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    aspect_ratio = 1,
    grid = true,
    framestyle = :box,
    )

    Colorbar(F[2,1:2], hm, vertical = false)
    return F
end

"""
2b) alternative heatmap: SSA vs NLSA
#take the H from the saved files
"""

function heatmap_comparison_clipped(F,ssa_h,nlsa_h,highclip,fs,IGBP_list,spotslist,variables_names)

    nlsa_h[nlsa_h .< 2] .= NaN
    ssa_h[ssa_h .< 2] .= NaN

    ax = Axis(F[1,1],
    #aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    #xticklabelalign = (:left,:top),
    xticklabelsize = fs,
    yticks = (1:length(variables_names),variables_names),
    #yticklabelalign = (:right,:center),
    yticklabelsize = fs,
    subtitle = "ssa"
    )



    hm = heatmap!(ax,ssa_h,
    colormap = cgrad(:heat, highclip, categorical = true),
    colorrange = (1,highclip),
    highclip = :red,
    lowclip = :white,
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    #aspect_ratio = 1,
    grid = true,
    framestyle = :box,

    )

    for i in 1:size(ssa_h)[1]
        str = ssa_h[i,:]
        map(enumerate(str)) do (j, c)
            if !isnan(c)
                text!(ax,i,j,text="$(Int(c))",
                align = (:center,:center),
                color = :black,fontsize=fs,
                )
            end
            #text!(ax,j,i,text=rich("$(lpad(string(c),padnr,"0"))",
            #text!(ax,i,j,text="$c",
            #color = :black,fontsize=fs)
        end
        #Label(f[1+i,1], rich(row_color_chars...),fontsize=18)
    end

    ax = Axis(F[1,2],
    #aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    #xticklabelalign = (:left,:top),
    xticklabelsize = fs,
    yticks = (1:length(variables_names),variables_names),
    #yticklabelalign = (:right,:center),
    yticklabelsize = fs,
    subtitle = "nlsa"
    )

    hm = heatmap!(ax,nlsa_h,
    colormap = cgrad(:heat, highclip, categorical = true),
    colorrange = (1,highclip),
    highclip = :red,
    lowclip = :white,
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    #aspect_ratio = 1,
    grid = true,
    framestyle = :box,
    )


    for i in 1:size(nlsa_h)[1]
        str = nlsa_h[i,:]
        map(enumerate(str)) do (j, c)
            if !isnan(c)
                text!(ax,i,j,text="$(Int(c))",
                align = (:center,:center),
                color = :black,fontsize=fs,
                )
            end
            #text!(ax,j,i,text=rich("$(lpad(string(c),padnr,"0"))",
            #text!(ax,i,j,text="$c",
            #color = :black,fontsize=fs)
        end
        #Label(f[1+i,1], rich(row_color_chars...),fontsize=18)
    end


    Colorbar(F[2,1:2], hm, vertical = false)
    return F
end

"""
3) lat long igbp
"""


function plot_errors(ax,yeartable,x,textstrings,color)



    y = mean(yeartable,dims=1)[:]
    ye = std(yeartable,dims=1)[:]

    scatter!(ax,x[:],y[:],linewidth=2,
    color=:black,marker=:x,markersize=10)

    errorbars!(ax,x,y,ye, color = color)

    for (i,xi) in enumerate(x)
        text!(
            ax, xi, y[i],
            text = textstrings[i], 
            align = (:center, :bottom),
            offset = (0, 2),
            #space = :relative,
            fontsize = 12
        )
    end
    return f
end

function plot_error_variables_loader(f,spots_list,vari)

    Fi = load("/net/home/lschulz/logs/KW_2_15/trends.jld2")
    ssa_trends = Fi["ssa_trends"]
    ssa_h = Fi["ssa_h"]
    nlsa_h = Fi["nlsa_h"]
    nlsa_trends = Fi["nlsa_trends"]


    # Calculate maximum year from the length of the time series
    maxyear = Int(floor(N/365))

    # Extract harmonic trends for ssa and nlsa modes for all selected years
    ssa_harm = [ssa_trends[spot,vari,i] for i in year_ind, spot in spots_list]
    nlsa_harm = [nlsa_trends[spot,vari,i] for i in year_ind, spot in spots_list]

    ssa_maxs = [maximum(ssa_harm[i,j]) for i in 1:maxyear, j in 1:length(spots_list)]
    nlsa_maxs = [maximum(nlsa_harm[i,j]) for i in 1:maxyear, j in 1:length(spots_list)]

    ssa_argmaxs = [argmax(ssa_harm[i,j]) for i in 1:maxyear, j in 1:length(spots_list)]
    nlsa_argmaxs = [argmax(nlsa_harm[i,j]) for i in 1:maxyear, j in 1:length(spots_list)]

    return ssa_maxs,nlsa_maxs,ssa_argmaxs,nlsa_argmaxs

end

function plot_error_wrapper_enf()
    inds = findall(x->x=="ENF",IGBP_list)
    spots_list = inds
    for vari = 1:16
        f = Figure(resolution=(900,1200))
        ssa_maxs,nlsa_maxs,ssa_argmaxs,nlsa_argmaxs = plot_error_variables_loader(f,spots_list,vari)


        ax = Axis(f[1,1],
        title = "amplitude lattitude dependency",
        xlabel = "lattitude",
        ylabel = "seasonal max",
        )


        stringlist = [spotslist[i].* "\n$(elevation[i]) m" for i in inds]

        plot_errors(ax,ssa_maxs,lattitude[spots_list],stringlist,color_ssa)
        plot_errors(ax,nlsa_maxs,lattitude[spots_list],stringlist,color_nlsa)

        ax = Axis(f[1,2],
        title = "phase lattitude dependency",
        xlabel = "lattitude",
        ylabel = "max position",
        )

        plot_errors(ax,ssa_argmaxs,lattitude[spots_list],stringlist,color_ssa)
        plot_errors(ax,nlsa_argmaxs,lattitude[spots_list],stringlist,color_nlsa)

        ax = Axis(f[2,1],
        title = "amplitude longitude dependency",
        xlabel = "longitude",
        ylabel = "seasonal max",
        )

        plot_errors(ax,ssa_maxs,longitude[spots_list],stringlist,color_ssa)
        plot_errors(ax,nlsa_maxs,longitude[spots_list],stringlist,color_nlsa)

        ax = Axis(f[2,2],
        title = "phase longitude dependency",
        xlabel = "longitude",
        ylabel = "max position",
        )

        plot_errors(ax,ssa_argmaxs,longitude[spots_list],stringlist,color_ssa)
        plot_errors(ax,nlsa_argmaxs,longitude[spots_list],stringlist,color_nlsa)
        
        ax = Axis(f[3,1],
        title = "amplitude elevation dependency",
        xlabel = "elevation /m",
        ylabel = "seasonal max",
        )
        stringlist = [spotslist[i].* "\n$(lattitude[i])\n$(longitude[i])" for i in inds]
        plot_errors(ax,ssa_maxs,elevation[spots_list],stringlist,color_ssa)
        plot_errors(ax,nlsa_maxs,elevation[spots_list],stringlist,color_nlsa)

        ax = Axis(f[3,2],
        title = "phase elevation dependency",
        xlabel = "elevation /m",
        ylabel = "max position",
        )

        plot_errors(f[3,2],ssa_argmaxs,elevation[spots_list],stringlist,color_ssa)
        plot_errors(f[3,2],nlsa_argmaxs,elevation[spots_list],stringlist,color_nlsa)

        save(dir*"locality_enf/$vari-$(variables_names[vari]).png",f)
    end
    return nothing
end

"""
3)b alternative lat long phase amplitude dependency
"""


function plot_errors(ax,yeartable,x,textstrings,color)

    years = startyear:startyear+Int(floor(N/365)-1)

    #y = maximum(yeartable,dims=3)[:]

    #scatter!(ax,years,yeartable,linewidth=2,
    #color=:black,marker=:x,markersize=10)

    series!(ax,years,yeartable',linewidth=2,
    solid_color=color,marker=:x,markersize=10)

    #errorbars!(ax,x,y,ye, color = color)

    
    for (i,xi) in enumerate(x)
        j = argmax(yeartable[:,i])
        text!(
            ax, years[j], yeartable[j,i],
            text = textstrings[i], 
            align = (:center, :bottom),
            offset = (0, -50),
            #space = :relative,
            fontsize = 12
        )
    end
    
    return ax
end

function plot_error_wrapper_enf()
    inds = findall(x->x=="ENF",IGBP_list)
    spots_list = inds
    for vari = 1:16
        f = Figure(resolution=(900,1200))

        ssa_maxs,nlsa_maxs,ssa_argmaxs,nlsa_argmaxs = plot_error_variables_loader(f,spots_list,vari)

        stringlist = [spotslist[i].* "\n$(lattitude[i])\n$(longitude[i])\n$(elevation[i]) m" for i in inds]

        year_offset = 2

        ax = Axis(f[1,1],
        title = "amplitude change ssa",
        xlabel = "time /a",
        ylabel = "seasonal max",
        limits=(startyear-year_offset,startyear+Int(floor(N/365)-1)+year_offset,0.5,2),
        )

        plot_errors(ax,ssa_maxs,lattitude[spots_list],stringlist,color_ssa)

        ax = Axis(f[1,2],
        title = "amplitude change nlsa",
        xlabel = "time /a",
        ylabel = "seasonal max",
        limits=(startyear-year_offset,startyear+Int(floor(N/365)-1)+year_offset,0.5,2),

        )

        plot_errors(ax,nlsa_maxs,lattitude[spots_list],stringlist,color_nlsa)


        ax = Axis(f[2,1],
        title = "phase change ssa",
        xlabel = "time /a",
        ylabel = "max position",
        limits=(startyear-year_offset,startyear+Int(floor(N/365)-1)+year_offset,120,220),

        )

        plot_errors(ax,ssa_argmaxs,lattitude[spots_list],stringlist,color_ssa)


        ax = Axis(f[2,2],
        title = "phase change nlsa",
        xlabel = "time /a",
        ylabel = "max position",
        limits=(startyear-year_offset,startyear+Int(floor(N/365)-1)+year_offset,120,220),

        )

        plot_errors(ax,nlsa_argmaxs,lattitude[spots_list],stringlist,color_nlsa)

        save(dir*"locality_enf/alternative_$vari-$(variables_names[vari]).png",f)
    end
    return nothing
end


"""
4) vegetation response plots
"""

function plot_veg_response(f,spots_list,vari_1,vari_2,stringlist,fontsize)

    Fi = load("/net/home/lschulz/logs/KW_2_15/trends.jld2")
    ssa_trends = Fi["ssa_trends"]
    ssa_h = Fi["ssa_h"]
    nlsa_h = Fi["nlsa_h"]
    nlsa_trends = Fi["nlsa_trends"]

    lags = -45:45
    xtiks = (-40:20:40,string.(-40:20:40))

    signal_1 = wholedata[:,spots_list,vari_1]'
    signal_2 = wholedata[:,spots_list,vari_2]'

    season_1_ssa = ssa_trends[spots_list,vari_1,:]
    season_2_ssa = ssa_trends[spots_list,vari_2,:]

    season_1_nlsa = nlsa_trends[spots_list,vari_1,:]
    season_2_nlsa = nlsa_trends[spots_list,vari_2,:]

    residual_1_ssa = signal_1 .- season_1_ssa
    residual_2_ssa = signal_2 .- season_2_ssa

    residual_1_nlsa = signal_1 .- season_1_nlsa
    residual_2_nlsa = signal_2 .- season_2_nlsa

    signal_cc = hcat([crosscor(signal_1[spot,:],signal_2[spot,:],lags) for spot in 1:length(spots_list)]...)
    season_ssa_cc = hcat([crosscor(season_1_ssa[spot,:],season_2_ssa[spot,:],lags) for spot in 1:length(spots_list)]...)
    season_nlsa_cc = hcat([crosscor(season_1_nlsa[spot,:],season_2_nlsa[spot,:],lags) for spot in 1:length(spots_list)]...)
    residual_ssa_cc = hcat([crosscor(residual_1_ssa[spot,:],residual_2_ssa[spot,:],lags) for spot in 1:length(spots_list)]...)
    residual_nlsa_cc = hcat([crosscor(residual_1_nlsa[spot,:],residual_2_nlsa[spot,:],lags) for spot in 1:length(spots_list)]...)

    season_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_1]] .= NaN
    season_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_2]] .= NaN
    season_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_1]] .= NaN
    season_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_2]] .= NaN
    residual_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_1]] .= NaN
    residual_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_2]] .= NaN
    residual_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_1]] .= NaN
    residual_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_2]] .= NaN


    ylim1 = 0.4
    gridwidth = 3

    ax1 = Axis(f[1,1],
    title = "signal",
    titlesize = fontsize,
    xlabel = "lag / d",
    ylabel = "cross correlation",
    limits=(lags[1],lags[end],ylim1,1),
    yticks = (ylim1:0.1:1,string.(ylim1:0.1:1)),
    xticks = xtiks,
    yticklabelsize = fontsize-4,
    ylabelsize = fontsize,
    xticklabelsize = fontsize-4,
    xlabelsize = fontsize,
    xgridwidth = gridwidth,
    ygridwidth = gridwidth,
    )


    ax2 = Axis(f[1,2],
    title = "season",
    titlesize = fontsize,
    xlabel = "lag / d",
    ylabel = "cross correlation",
    limits=(lags[1],lags[end],ylim1,1),
    yticks = (ylim1:0.1:1,string.(ylim1:0.1:1)),
    xticks = xtiks,
    yticklabelsize = fontsize-4,
    ylabelsize = fontsize,
    xticklabelsize = fontsize-4,
    xlabelsize = fontsize,
    xgridwidth = gridwidth,
    ygridwidth = gridwidth,
    )

    ax3 = Axis(f[1,3],
    title = "residual",
    titlesize = fontsize,
    xlabel = "lag / d",
    ylabel = "cross correlation",
    limits=(lags[1],lags[end],ylim1,1),
    yticks = (ylim1:0.1:1,string.(ylim1:0.1:1)),
    xticks = xtiks,
    yticklabelsize = fontsize-4,
    ylabelsize = fontsize,
    xticklabelsize = fontsize-4,
    xlabelsize = fontsize,
    xgridwidth = gridwidth,
    ygridwidth = gridwidth,
    )

    series!(ax1, lags, signal_cc', solid_color = color_signal)
    series!(ax2, lags, season_ssa_cc', solid_color = color_ssa)
    series!(ax3, lags, residual_ssa_cc', solid_color = color_ssa)
    series!(ax2, lags, season_nlsa_cc', solid_color = color_nlsa)
    series!(ax3, lags, residual_nlsa_cc', solid_color = color_nlsa)

    hideydecorations!(ax2,grid = false)
    hideydecorations!(ax3,grid = false)

    return f
end

function vegetation_response_enf()
    inds = findall(x->x=="ENF",IGBP_list)
    fontsize = 26

    for (i,j) in [[1,16],[1,12],[1,13],[1,14],[2,16],[2,12],[2,13],[2,14],
        [3,16],[3,12],[3,13],[3,14],[4,16],[4,12],[4,13],[4,14]]
    f = plot_veg_response(Figure(resolution=(800,400)),inds,i,j,[""],fontsize)
    v1 = variables_names[i]
    v2 = variables_names[j]
    save(dir*"response/veg_response_enf_$(v1)_$v2.png",f)
    end
end

"""
mask for selection 
"""
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

#1. 4/a #2. 7/a
# masks: igbp, vari, 
# 1.     yes    no
# 2.     no     yes
# 3.     yes    yes
# 4.     no     no

function apply_masks_do_harmonics()
    ssa_h,nlsa_h,ssa_trends,nlsa_trends = create_trend_tensors(N,W)

    for mask1= [0,1], mask2 = [0,1]

        variables_names = ["GPP_DAY_1","GPP_DAY_2","GPP_NIGHT_1","GPP_NIGHT_2",
        "RECO_DAY_1","RECO_DAY_2","RECO_NIGHT_1","RECO_NIGHT_2",
        "NEE","NEE_NIGHT","NEE_DAY",
        "TA","SW_IN","LW_IN","SWC","TS"]

        #medium length
        spotslist = [
            "BE-Lon",
            "BE-Vie",
            "CH-Dav",
            "CH-Fru",
            "CH-Lae",
            "CH-Oe2",
            "DE-Geb",
            "DE-Gri",
            "DE-Hai",
            "DE-Tha",
            "DK-Sor",
            "ES-LJu",
            "FI-Hyy",
            "FR-Aur",
            "FR-Lam",
            "IT-BCi",
            "IT-Lav",
            "RU-Fyo",
        ]

        IGBP_list = [
            "CRO",
            "MF",
            "ENF",
            "GRA",
            "MF",
            "CRO",
            "CRO",
            "GRA",
            "DBF",
            "ENF",
            "DBF",
            "OSH",
            "ENF",
            "CRO",
            "CRO",
            "CRO",
            "ENF",
            "ENF",
        ]

        vari_list = 1:length(variables_names)

        highclip = 9
        fs = 8

        if mask1 == true
            forest_list,grass_list = mask_IGBP(IGBP_list)
            spots_list = forest_list
            spotslist = spotslist[spots_list]
            IGBP_list = IGBP_list[forest_list]
        else
            spots_list = 1:length(spotslist)
        end
        if mask2 == true
            vari_list,variables_names = mask_vari(variables_names)
        else
            vari_list = 1:length(variables_names)
        end

        F = Figure()
        println(size(ssa_h[spots_list,vari_list]))
        heatmap_comparison_clipped(F,ssa_h[spots_list,vari_list],nlsa_h[spots_list,vari_list],highclip,fs,IGBP_list,spotslist,variables_names)
        save(dir*"raw_IGBP_$(mask1)_vari_$(mask2).png",F)
    end



end


#iterate combined plot over the spots and variables
function plot_existing_data_masked()


    variables_names = ["GPP_DAY_1","GPP_DAY_2","GPP_NIGHT_1","GPP_NIGHT_2",
    "RECO_DAY_1","RECO_DAY_2","RECO_NIGHT_1","RECO_NIGHT_2",
    "NEE","NEE_NIGHT","NEE_DAY",
    "TA","SW_IN","LW_IN","SWC","TS"]

    #medium length
    spotslist = [
        "BE-Lon",
        "BE-Vie",
        "CH-Dav",
        "CH-Fru",
        "CH-Lae",
        "CH-Oe2",
        "DE-Geb",
        "DE-Gri",
        "DE-Hai",
        "DE-Tha",
        "DK-Sor",
        "ES-LJu",
        "FI-Hyy",
        "FR-Aur",
        "FR-Lam",
        "IT-BCi",
        "IT-Lav",
        "RU-Fyo",
    ]

    IGBP_list = [
        "CRO",
        "MF",
        "ENF",
        "GRA",
        "MF",
        "CRO",
        "CRO",
        "GRA",
        "DBF",
        "ENF",
        "DBF",
        "OSH",
        "ENF",
        "CRO",
        "CRO",
        "CRO",
        "ENF",
        "ENF",
    ]


    forest_list,grass_list = mask_IGBP(IGBP_list)
    spots_list = forest_list
    #spotslist = spotslist[spots_list]
    #IGBP_list = IGBP_list[forest_list]
    vari_list,variables_names_2 = mask_vari(variables_names)

    for vari = vari_list, spot = spots_list
        try
            savedirname = dir*"runs/spot$spot/"
            create_combined_plot(spot,savedirname,W,vari,N)
            println("YES spot $spot variable $vari")
        catch
            println("NO spot $spot variable $vari")
        end
    end
end

"""
heatmap for paper
    ssa_h,nlsa_h,ssa_trends,nlsa_trends = create_trend_tensors(N,W)
"""

function heatmap_paper(ssa_h,nlsa_h)

    variables_names = ["GPP_DAY_1","GPP_DAY_2","GPP_NIGHT_1","GPP_NIGHT_2",
    "RECO_DAY_1","RECO_DAY_2","RECO_NIGHT_1","RECO_NIGHT_2",
    "NEE","NEE_NIGHT","NEE_DAY",
    "TA","SW_IN","LW_IN","SWC","TS"]

    #medium length
    spotslist = [
        "BE-Lon",
        "BE-Vie",
        "CH-Dav",
        "CH-Fru",
        "CH-Lae",
        "CH-Oe2",
        "DE-Geb",
        "DE-Gri",
        "DE-Hai",
        "DE-Tha",
        "DK-Sor",
        "ES-LJu",
        "FI-Hyy",
        "FR-Aur",
        "FR-Lam",
        "IT-BCi",
        "IT-Lav",
        "RU-Fyo",
    ]

    IGBP_list = [
        "CRO",
        "MF",
        "ENF",
        "GRA",
        "MF",
        "CRO",
        "CRO",
        "GRA",
        "DBF",
        "ENF",
        "DBF",
        "OSH",
        "ENF",
        "CRO",
        "CRO",
        "CRO",
        "ENF",
        "ENF",
    ]
    
    forest_list,grass_list = mask_IGBP(IGBP_list)
    spots_list = forest_list
    spotslist = spotslist[spots_list]
    IGBP_list = IGBP_list[forest_list]
    vari_list,variables_names = mask_vari(variables_names)

    F = Figure(resolution = (800,500))
    
    highclip = 8
    fs = 26

    nlsa_h[nlsa_h .< 2] .= NaN
    ssa_h[ssa_h .< 2] .= NaN

    ax = Axis(F[1,1],
    #aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    xlabel="site",
    #ylabel = "variable",
    xlabelsize = fs,
    xticklabelsize = fs-4,
    yticks = (1:length(variables_names),variables_names),
    yticklabelsize = fs-4,
    ylabelsize = fs,
    title = "SSA",
    titlesize = fs,
    )



    hm = heatmap!(ax,ssa_h[spots_list,vari_list],
    colormap = cgrad(:heat, highclip, categorical = true),
    colorrange = (1,highclip),
    highclip = :red,
    lowclip = :white,
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    #aspect_ratio = 1,
    grid = true,
    framestyle = :box,

    )

    for i in 1:size(ssa_h[spots_list,vari_list])[1]
        str = ssa_h[spots_list,vari_list][i,:]
        map(enumerate(str)) do (j, c)
            if !isnan(c)
                text!(ax,i,j,text="$(Int(c))",
                align = (:center,:center),
                color = :black,fontsize=fs,
                )
            end
            #text!(ax,j,i,text=rich("$(lpad(string(c),padnr,"0"))",
            #text!(ax,i,j,text="$c",
            #color = :black,fontsize=fs)
        end
        #Label(f[1+i,1], rich(row_color_chars...),fontsize=18)
    end

    ax = Axis(F[1,2],
    #aspect = AxisAspect(1),
    xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    xlabel="site",
    #ylabel = "variable",
    xlabelsize = fs,
    xticklabelsize = fs-4,
    yticks = (1:length(variables_names),variables_names),
    yticklabelsize = fs-4,
    ylabelsize = fs,
    title = "NLSA",
    titlesize = fs,
    )

    hm = heatmap!(ax,nlsa_h[spots_list,vari_list],
    colormap = cgrad(:heat, highclip, categorical = true),
    colorrange = (1,highclip),
    highclip = :red,
    lowclip = :white,
    #colorrange = (minimum(valuetable),maximum(valuetable)),
    #aspect_ratio = 1,
    grid = true,
    framestyle = :box,
    )


    for i in 1:size(nlsa_h[spots_list,vari_list])[1]
        str = nlsa_h[spots_list,vari_list][i,:]
        map(enumerate(str)) do (j, c)
            if !isnan(c)
                text!(ax,i,j,text="$(Int(c))",
                align = (:center,:center),
                color = :black,fontsize=fs,
                )
            end
            #text!(ax,j,i,text=rich("$(lpad(string(c),padnr,"0"))",
            #text!(ax,i,j,text="$c",
            #color = :black,fontsize=fs)
        end
        #Label(f[1+i,1], rich(row_color_chars...),fontsize=18)
    end

    hideydecorations!(ax)

    Colorbar(F[2,1:2], hm, vertical = false,label = "Harmonics for raw signal",labelsize = fs)

    save(dir*"heatmap_harmonics_raw.png",F)
end

"""
mode figure for paper
"""

function mode_figure_paper(p)

    p = rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    modenumber = 6
    smallfs = 16
    fs = 26

    F = Figure(resolution=(800,800))

    #figure

    ax_time = Axis(F[1,1:4],xticks=Int.(floor.(years[1]:3:years[end])),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time (a)",
    ylabel=varname,
    xlabelsize = fs,
    xticklabelsize = fs-4,
    yticklabelsize = fs-4,
    ylabelsize = fs,)

    hideydecorations!(ax_time)

    ax_spec = Axis(F[2,1:4],yscale=log10,
    xlabel="frequency (/a)",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs-4,
    yticklabelsize = fs-4,
    ylabelsize = fs,)

    hideydecorations!(ax_spec)

    ax_modes = Axis(F[3:4,1:2],yticksvisible = false,
    yticklabelsvisible = false,xticks=1:Int(floor(W/365.25)),
    xlabel="time (a)",
    ylabel="individual modes",
    xlabelsize = fs,
    xticklabelsize = fs-4,
    yticklabelsize = fs-4,
    ylabelsize = fs,)

    hideydecorations!(ax_modes)

    ax_freq = Axis(F[3:4,3:4],limits=(1/10, 7,-1,modenumber*2+3),#,yscale=log10
    xticks=1:7,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency (/a)",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs-4,
    yticklabelsize = fs-4,
    ylabelsize = fs,)

    hideydecorations!(ax_freq)

    #ax_labels = Axis(F[3:4,5])
    #hidedecorations!(ax_labels)

    #plotting
    scatter!(ax_time,years,signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_time,years,signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_time,years,ssa_rec,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_time,years,nlsa_rec,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain_N,spec_signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_spec,freq_domain_N,spec_ssa_rc,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain_N,spec_nlsa_rc,color=color_nlsa,linewidth=lw,label="NLSA")

    axislegend(ax_spec)

    text!(
        ax_time, 0, 1,
        text = "signal and reconstructions", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_spec, 0, 1,
        text = "spectra of signal and reconstructions", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_time, 0.95, 0.95,
        text = "$(spotslist[spot])\n$(igbpclass)\n$(varname)", 
        align = (:right, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    #modes 

    color_fit = "grey"
    for i = 1:modenumber

        mode = ssa_Eof[:,i]
        Four = spec_ssa_eof[:,i]
        #Four ./= maximum(abs.(Four)) 

        years_modes = (1:length(mode)) ./ 365.25


        lines!(ax_modes,years_modes,mode .+ (i*0.07),
        color=color_ssa)
        lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_ssa[i]) .+ i,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i, #[freqstart:end]
        color=color_ssa)

        mode = nlsa_Eof[:,i] 
        Four = spec_nlsa_eof[:,i]


        lines!(ax_modes,years_modes,mode .+ ((modenumber +1+i)*0.08),
        color=color_nlsa)
        lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_nlsa[i]) .+ i .+ modenumber .+1,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i .+ modenumber .+1, #[freqstart:end]
        color=color_nlsa)

    end

    #height and freq
    boxstring = "f \t height \n"#"peaks \n\n f \t height \n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_nlsa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_nlsa[k][2],digits=1))
        boxstring *= "\n"
    end
    boxstring*= "\n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_ssa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_ssa[k][2],digits=1))
        boxstring *= "\n"
    end

    """
    text!(
        ax_labels, 0, 1,
        text = boxstring, 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = smallfs
    )
    """ 

    text!(
        ax_modes, 0, 1,
        text = "mode shapes", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_freq, 0, 1,
        text = "mode spectra", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )


    #plot_seasonality(F[2,1],p)
    #plot_annuality(F[2,2],p)

    save(dir*"modes_raw.png",F)

end

"""
introduction figure paper
"""

function intro_figure_paper(p,pic_name,varname_resolved)

    #set_theme!(fonts = (; regular = "Computer Modern", bold = "Computer Modern"))

    pic_name = "modes/"*pic_name
    
    p = rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    modenumber = 8
    modenumber_k = 48
    xstep_k = 16
    k_xticks = (xstep_k:xstep_k:modenumber_k,string.(xstep_k:xstep_k:modenumber_k))
    spec_xticks = (1:5,string.(1:5))
    smallfs = 16
    fs = 23
    speclimits = (freq_domain_N[1],6,10^-4,1)
    klimits = (0,modenumber_k+1,10^-4,1)
    modeslimits = (0,2556/365.25,-0.04,(2*modenumber+2.5)*0.08)
    modeyticks = (0.08:0.08*2:(2*modenumber+1)*0.08,string.(vcat(1:2:modenumber,"",1:2:modenumber)))
    freqlimits = (1/10, 7,-0.2,modenumber*2+2.5)
    gapsize = 8
    lw = 4
    ms = 12
    scatterstep = 10
    color_signal = "grey68"

    F = Figure(resolution=(800,800))

    #figure

    ax_time = Axis(F[1,1:3],
    xticks=Int.(floor.(years[1]:3:years[end])),
    limits = (years[1],years[end],minimum(signal)*1.1,maximum(signal)*1.1),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time [a]",
    ylabel=varname_resolved,
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    #hideydecorations!(ax_time,label = false)

    ax_spec = Axis(F[2,1:3],yscale=log10,
    limits = speclimits,
    xticks = spec_xticks,
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    #hideydecorations!(ax_spec,label = false)

    ax_k = Axis(F[2,4],yscale=log10,
    limits = klimits,
    xlabel="mode dim",
    xlabelsize = fs,
    xticks = k_xticks,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_k,grid = false)

    ax_modes = Axis(F[3:4,1:2],#yticksvisible = false,
    #yticklabelsvisible = false,
    xticks=0:Int(floor(W/365.25)),
    yticks = modeyticks,
    limits = modeslimits,
    xlabel="time [a]",
    ylabel="individual modes",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    #hideydecorations!(ax_modes,label=false)

    ax_freq = Axis(F[3:4,3:4],
    limits=freqlimits,#,yscale=log10
    xticks=1:7,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_freq)

    #ax_labels = Axis(F[3:4,5])
    #hidedecorations!(ax_labels)

    #plotting
    #scatter!(ax_time,years[1:scatterstep:end],signal[1:scatterstep:end],color=color_signal,markersize = ms,marker=:x)
    signal_l = lines!(ax_time,years,signal,color=color_signal,linewidth=lw,label = "signal")
    ssa_l = lines!(ax_time,years,ssa_rec,color=color_ssa,linewidth=lw,label="SSA")
    nlsa_l = lines!(ax_time,years,nlsa_rec,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain_N,spec_signal,color=color_signal,linewidth=lw,label = "signal")
    signal_s = scatter!(ax_spec,freq_domain_N,spec_signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_spec,freq_domain_N,spec_ssa_rc,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain_N,spec_nlsa_rc,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_k,1:modenumber_k,ssa_cap_var[1:modenumber_k],color=color_ssa,linewidth=lw)
    lines!(ax_k,1:modenumber_k,nlsa_cap_var[1:modenumber_k],color=color_nlsa,linewidth=lw)




    #modes 

    color_fit = "grey"
    for i = 1:modenumber

        mode = ssa_Eof[:,i]
        Four = spec_ssa_eof[:,i]
        #Four ./= maximum(abs.(Four)) 

        years_modes = (1:length(mode)) ./ 365.25


        lines!(ax_modes,years_modes,mode .+ (i*0.08),
        color=color_ssa)
        lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_ssa[i]) .+ i,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i, #[freqstart:end]
        color=color_ssa)

        mode = nlsa_Eof[:,i] 
        Four = spec_nlsa_eof[:,i]


        lines!(ax_modes,years_modes,mode .+ ((modenumber +1+i)*0.08),
        color=color_nlsa)
        lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_nlsa[i]) .+ i .+ modenumber .+1,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i .+ modenumber .+1, #[freqstart:end]
        color=color_nlsa)

    end

    #height and freq
    boxstring = "f \t height \n"#"peaks \n\n f \t height \n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_nlsa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_nlsa[k][2],digits=1))
        boxstring *= "\n"
    end
    boxstring*= "\n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_ssa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_ssa[k][2],digits=1))
        boxstring *= "\n"
    end

    text!(
        ax_time, 0, 1,
        #text = "signal and reconstructions", 
        text = "(a)",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_spec, 0, 1,
        #text = "signal and reconstructions", 
        text = "(b)",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_k, 0.7, 1,
        #text = "signal and reconstructions", 
        text = "(c)",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_modes, 0, 1,
        text = "(d)", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax_freq, 0, 1,
        text = "(e)", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    linkyaxes!(ax_spec, ax_k)





    Legend(F[1,4],
    #[lin, sca, [lin, sca], sca2],
    #["a line", "some dots", "both together", "rect markers"])
    [[signal_l,signal_s],ssa_l,nlsa_l],
    ["signal","SSA","NLSA"],
    "$(spotslist[spot])\n$(igbpclass)",
    framevisible = false,
    labelsize = fs,
    titlesize=fs+4,
    #titlehalign = :left,
    #titlevalign = :bottom,
    #titleposition = :left,
    )
    colgap!(F.layout, 1, gapsize) #colgap!(fig.layout, 1, Relative(0.15))
    colgap!(F.layout, 2, gapsize)
    colgap!(F.layout, 3, gapsize)
    rowgap!(F.layout, 1, gapsize)
    rowgap!(F.layout, 2, gapsize)
    rowgap!(F.layout, 3, gapsize)




    save(dir*pic_name*".png",F)

end

function create_paper_overview_plots()
    """raw"""
    outdir="/net/scratch/lschulz/fluxfullset_midwithnee/"

    varname_r = "GPP"
    spot = 9
    vari = 1
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"raw_intro_$(spotname)_$(varname)",varname_r)

    varname_r = "TS"
    spot = 5
    vari = 16
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"raw_resolved_$(spotname)_$(varname)",varname_r)

    varname_r = "RECO"
    spot = 5
    vari = 7
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"raw_unresolved_$(spotname)_$(varname)",varname_r)

    varname_r = "TS"
    spot = 17
    vari = 16
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"raw_measurement_$(spotname)_$(varname)",varname_r)

    """filter"""
    outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass7/"

    varname_r = "GPP"
    spot = 9
    vari = 1
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"filter_intro_$(spotname)_$(varname)",varname_r)

    varname_r = "TS"
    spot = 5
    vari = 16
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"filter_resolved_$(spotname)_$(varname)",varname_r)

    varname_r = "RECO"
    spot = 5
    vari = 7
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"filter_unresolved_$(spotname)_$(varname)",varname_r)

    varname_r = "TS"
    spot = 17
    vari = 16
    spotname = spotslist[spot]
    varname = variables_names[vari]
    p = local_parameters(W,vari,spot,N,startyear)
    intro_figure_paper(p,"filter_measurement_$(spotname)_$(varname)",varname_r)
end

"""
crosscorrelation figure paper
"""

function plot_veg_response(vari_1,vari_2,savename)

    spot = 18
    spots_list = [17,18,3,10,13]
    spots_list = [3]
    #vari_1 = 1
    #vari_2 = 12
    stringlist = "test"
    fontsize = 23
    #savename = "t/$vari_2"
    fs = 18

    ylim1 = 0.4
    gridwidth = 3
    fontsize_plot = 10

    lags = -60:60
    xtiks = (-60:20:60,string.(-60:20:60))


    f = Figure()

    Fi = load("/net/home/lschulz/logs/KW_2_15/trends.jld2")
    ssa_trends = Fi["ssa_trends"]
    ssa_h = Fi["ssa_h"]
    nlsa_h = Fi["nlsa_h"]
    nlsa_trends = Fi["nlsa_trends"]



    signal_1 = wholedata[:,spots_list,vari_1]'
    signal_2 = wholedata[:,spots_list,vari_2]'

    season_1_ssa = ssa_trends[spots_list,vari_1,:]
    season_2_ssa = ssa_trends[spots_list,vari_2,:]

    season_1_nlsa = nlsa_trends[spots_list,vari_1,:]
    season_2_nlsa = nlsa_trends[spots_list,vari_2,:]

    residual_1_ssa = signal_1 .- season_1_ssa
    residual_2_ssa = signal_2 .- season_2_ssa

    residual_1_nlsa = signal_1 .- season_1_nlsa
    residual_2_nlsa = signal_2 .- season_2_nlsa

    signal_cc = hcat([crosscor(signal_1[spot,:],signal_2[spot,:],lags) for spot in 1:length(spots_list)]...)
    season_ssa_cc = hcat([crosscor(season_1_ssa[spot,:],season_2_ssa[spot,:],lags) for spot in 1:length(spots_list)]...)
    season_nlsa_cc = hcat([crosscor(season_1_nlsa[spot,:],season_2_nlsa[spot,:],lags) for spot in 1:length(spots_list)]...)
    residual_ssa_cc = hcat([crosscor(residual_1_ssa[spot,:],residual_2_ssa[spot,:],lags) for spot in 1:length(spots_list)]...)
    residual_nlsa_cc = hcat([crosscor(residual_1_nlsa[spot,:],residual_2_nlsa[spot,:],lags) for spot in 1:length(spots_list)]...)

    season_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_1]] .= NaN
    season_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_2]] .= NaN
    season_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_1]] .= NaN
    season_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_2]] .= NaN
    residual_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_1]] .= NaN
    residual_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_2]] .= NaN
    residual_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_1]] .= NaN
    residual_nlsa_cc[:,(nlsa_h .< 2)[spots_list,vari_2]] .= NaN


    ax1 = Axis(f[1,1],
    #title = "signal",
    titlesize = fontsize,
    xlabel = "lag / d",
    ylabel = "cross corr",
    limits=(lags[1],lags[end],ylim1,1.05),
    yticks = (ylim1:0.2:1,string.(ylim1:0.2:1)),
    xticks = xtiks,
    yticklabelsize = fontsize,
    ylabelsize = fontsize,
    xticklabelsize = fontsize,
    xlabelsize = fontsize,
    xgridwidth = gridwidth,
    ygridwidth = gridwidth,
    )


    series!(ax1, lags, signal_cc', solid_color = color_signal)
    series!(ax1, lags, season_ssa_cc', solid_color = color_ssa)
    series!(ax1, lags, residual_nlsa_cc', solid_color = color_nlsa)
    series!(ax1, lags, residual_ssa_cc', solid_color = color_ssa)
    series!(ax1, lags, season_nlsa_cc', solid_color = color_nlsa)

    amax_signal = argmax(signal_cc)
    amax_ssa = argmax(season_ssa_cc)
    amax_ssa_residual = argmax(residual_ssa_cc)
    amax_nlsa = argmax(season_nlsa_cc)
    amax_nlsa_residual = argmax(residual_nlsa_cc)

    scatter!(ax1,lags[amax_ssa],season_ssa_cc[amax_ssa],color = color_ssa,markersize = 10,marker=:x)
    scatter!(ax1,lags[amax_ssa_residual],residual_ssa_cc[amax_ssa_residual],color = color_ssa,markersize = 10,marker=:x)
    scatter!(ax1,lags[amax_signal],signal_cc[amax_signal],color = color_signal,markersize = 10,marker=:x)
    scatter!(ax1,lags[amax_nlsa],season_nlsa_cc[amax_nlsa],color = color_nlsa,markersize = 10,marker=:x)
    scatter!(ax1,lags[amax_nlsa_residual],residual_nlsa_cc[amax_nlsa_residual],color = color_nlsa,markersize = 10,marker=:x)


    #text!(ax1,lags[amax_signal],signal_cc[amax_signal],text="signal",align = (:center,:bottom),color = color_signal,fontsize=fontsize)
    #text!(ax1,lags[amax_ssa],season_ssa_cc[amax_ssa],text="SSA",align = (:center,:bottom),color = color_ssa,fontsize=fontsize)
    #text!(ax1,lags[amax_ssa_residual],residual_ssa_cc[amax_ssa_residual],text="SSA residual",align = (:center,:bottom),color = color_ssa,fontsize=fontsize)
    #text!(ax1,lags[amax_nlsa],season_nlsa_cc[amax_nlsa],text="NLSA",align = (:center,:bottom),color = color_nlsa,fontsize=fontsize)

    text!(
        ax1, 0.15, 0.75,
        #text = "signal and reconstructions", 
        text = "cycle",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax1, 0.4, 0.55,
        #text = "signal and reconstructions", 
        text = "signal",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )

    text!(
        ax1, 0.6, 0.3,
        #text = "signal and reconstructions", 
        text = "residual",
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fs
    )


    #series!(ax1, lags, season_nlsa_cc', solid_color = color_nlsa)
    #series!(ax1, lags, residual_nlsa_cc', solid_color = color_nlsa)

    #hideydecorations!(ax2,grid = false)
    #hideydecorations!(ax3,grid = false)

    ax2 = Axis(f[1,2:3],
    #title = "signal",
    titlesize = fontsize,
    xlabel = "lag / d",
    ylabel = "response value",
    #limits=(lags[1],lags[end],ylim1,1),
    #yticks = (ylim1:0.2:1,string.(ylim1:0.2:1)),
    #xticks = xtiks,
    yticklabelsize = fontsize,
    ylabelsize = fontsize,
    xticklabelsize = fontsize,
    xlabelsize = fontsize,
    xgridwidth = gridwidth,
    ygridwidth = gridwidth,
    )



    #linkyaxes!(ax1, ax2)
    hideydecorations!(ax2,ticks = false,ticklabels=false,grid = false)

    li_x = [0]
    li_y = []
    for spot = mask_IGBP(IGBP_list)[1]
        if nlsa_h[spot,vari_1] > 2 && nlsa_h[spot,vari_2] > 2

            text = spotslist[spot]*" "*IGBP_list[spot]
            season_1_ssa = ssa_trends[spot,vari_1,:]
            season_2_ssa = ssa_trends[spot,vari_2,:]

                
            season_1_nlsa = nlsa_trends[spot,vari_1,:]
            season_2_nlsa = nlsa_trends[spot,vari_2,:]

            season_ssa_cc = crosscor(season_1_ssa,season_2_ssa,lags)
            season_nlsa_cc = crosscor(season_1_nlsa,season_2_nlsa,lags)

            #season_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_1]] .= NaN
            #season_ssa_cc[:,(ssa_h .< 2)[spots_list,vari_2]] .= NaN

            amax_ssa = argmax(season_ssa_cc)
            amax_nlsa = argmax(season_nlsa_cc)

            x = lattitude[spot]
            value_ssa = season_ssa_cc[amax_ssa]
            value_nlsa = season_nlsa_cc[amax_nlsa]

            lag_ssa = lags[amax_ssa]
            lag_nlsa = lags[amax_nlsa]

            for v in [lag_ssa,lag_nlsa]
                if !isnan(v) && v != lags[1]
                    li_x = append!(li_x,v)
                end
            end

            for v in [value_ssa,value_nlsa]
                if !isnan(v)
                    li_y = append!(li_y,v)
                end
            end

            scatter!(ax2,lag_ssa,value_ssa,color = color_ssa,markersize = 10,marker=:x)
            scatter!(ax2,lag_nlsa,value_nlsa,color = color_nlsa,markersize = 10,marker=:x)

            text!(ax2,lag_ssa,value_ssa,text=text,
            align = (:center,:bottom),color = color_ssa,fontsize=fontsize_plot)

            try
                text!(ax2,lag_nlsa,value_nlsa,text=text,
                align = (:center,:bottom),color = color_nlsa,fontsize=fontsize_plot)
            catch
                println("not isnan $spot")
            end
        end
    end

    println(li_x)
    println(li_y)

    x1 = minimum(li_x)
    x2 = maximum(li_x)
    y1 = minimum(li_y)
    y2 = maximum(li_y)

    #println(x1)
    #println(x2)
    #println(y1)
    #println(y2)

    xpoints = [x1,x2,x2,x1,x1]
    ypoints = [y1,y1,y2,y2,y1]

    lines!(ax2,xpoints,ypoints,color=:black,linewidth=2)
    lines!(ax1,xpoints,ypoints,color=:black,linewidth=2)


    save(dir*savename*".png",f)
end

"""
plot_veg_response(13,1,"cc_SW-GPP")
plot_veg_response(13,5,"cc_SW-RECO")
plot_veg_response(13,9,"cc_SW-NEE")
plot_veg_response(13,12,"cc_SW-TA")
plot_veg_response(13,16,"cc_SW-TS")
plot_veg_response(12,1,"cc_TA-GPP")
plot_veg_response(16,1,"cc_TS-GPP")
"""

"""
modes figures
FONTS
set_theme!(fonts=(
    regular="Latin Modern Roman",
    bold = "Latin Modern Roman Bold",
    italic = "Latin Modern Roman Italic",
    bold_italic = "Latin Modern Roman Bold Italic",))

save(dir*"test.png",mode_figure(Figure(resolution=(800,1000)),p,"GPP"))
"""

function mode_figure(F,p,varname_resolved)

    p = rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    modenumber = 12
    modenumber_k = 48
    xstep_k = 16
    k_xticks = (xstep_k:xstep_k:modenumber_k,string.(xstep_k:xstep_k:modenumber_k))
    spec_xticks = (1:6,string.(1:6))
    smallfs = 16
    fs = 23
    speclimits = (freq_domain_N[1],6.05,10^-4,1)
    klimits = (0,modenumber_k+1,10^-4,1)
    modeslimits = (0,2556/365.25,-0.04,(2*modenumber+2.5)*0.08)
    modeyticks = (vcat(0.08:0.08*2:(modenumber)*0.08,(modenumber+2)*0.08:0.08*2:(2*modenumber+1)*0.08),string.(vcat(1:2:modenumber,1:2:modenumber)))
    freqlimits = (1/10, 7,-0.2,modenumber*2+2.5)
    gapsize = 8
    lw = 2
    ms = 12
    scatterstep = 10
    color_signal = "grey68"

    

    #figure

    ax_time = Axis(F[1:2,1:4],
    xticks=Int.(floor.(years[1]:3:years[end])),
    limits = (years[1],years[end],minimum(signal)*1.1,maximum(signal)*1.1),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time [a]",
    ylabel=varname_resolved,
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_time,ticks = false,ticklabels=false,grid=false,label=false)

    ax_spec = Axis(F[3:4,1:3],yscale=log10,
    limits = speclimits,
    xticks = spec_xticks,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(4),
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_spec,ticks = false,ticklabels=false,grid=false)

    ax_k = Axis(F[3:4,4],yscale=log10,
    limits = klimits,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(2),
    xlabel="dim",
    xlabelsize = fs,
    xticks = k_xticks,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_k,grid=false)


    ax_modes = Axis(F[5:9,1:2],#yticksvisible = false,
    #yticklabelsvisible = false,
    xticks=0:Int(floor(W/365.25)),
    yticks = modeyticks,
    limits = modeslimits,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(2),
    xlabel="time [a]",
    ylabel="individual modes",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_modes,ticks = false,ticklabels=false,grid=false)

    ax_freq = Axis(F[5:9,3:4],
    limits=freqlimits,#,yscale=log10
    xticks=1:7,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_freq)

    #ax_labels = Axis(F[3:4,5])
    #hidedecorations!(ax_labels)

    #plotting
    #scatter!(ax_time,years[1:scatterstep:end],signal[1:scatterstep:end],color=color_signal,markersize = ms,marker=:x)
    signal_l = lines!(ax_time,years,signal,color=color_signal,linewidth=lw,label = "signal")
    ssa_l = lines!(ax_time,years,ssa_trend_harm,color=color_ssa,linewidth=lw,label="SSA")
    nlsa_l = lines!(ax_time,years,nlsa_trend_harm,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain_N,spec_signal,color=color_signal,linewidth=lw,label = "signal")
    signal_s = scatter!(ax_spec,freq_domain_N,spec_signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_spec,freq_domain_N,spec_ssa_rc,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain_N,spec_nlsa_rc,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_k,1:modenumber_k,ssa_cap_var[1:modenumber_k],color=color_ssa,linewidth=lw)
    lines!(ax_k,1:modenumber_k,nlsa_cap_var[1:modenumber_k],color=color_nlsa,linewidth=lw)




    #modes 

    color_fit = "grey"
    for i = 1:modenumber

        mode = ssa_Eof[:,i]
        Four = spec_ssa_eof[:,i]
        #Four ./= maximum(abs.(Four)) 

        years_modes = (1:length(mode)) ./ 365.25


        lines!(ax_modes,years_modes,mode .+ (i*0.08),
        color=color_ssa)
        #lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_ssa[i]) .+ i,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i, #[freqstart:end]
        color=color_ssa)

        mode = nlsa_Eof[:,i] 
        Four = spec_nlsa_eof[:,i]


        lines!(ax_modes,years_modes,mode .+ ((modenumber +1+i)*0.08),
        color=color_nlsa)
        #lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_nlsa[i]) .+ i .+ modenumber .+1,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i .+ modenumber .+1, #[freqstart:end]
        color=color_nlsa)

    end

    #height and freq
    boxstring = "f \t height \n"#"peaks \n\n f \t height \n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_nlsa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_nlsa[k][2],digits=1))
        boxstring *= "\n"
    end
    boxstring*= "\n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_ssa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_ssa[k][2],digits=1))
        boxstring *= "\n"
    end
    linkyaxes!(ax_spec, ax_k)

    #colgap!(F.layout, 1, gapsize) #colgap!(fig.layout, 1, Relative(0.15))
    #colgap!(F.layout, 2, gapsize)
    #colgap!(F.layout, 3, gapsize)
    #rowgap!(F.layout, 1, gapsize)
    #rowgap!(F.layout, 2, gapsize)
    #rowgap!(F.layout, 3, gapsize)

    return F
end

function large_mode_figure(savedirname)
    F = Figure(resolution=(2400,1000))

    ga = F[1:12, 1:8] = GridLayout()
    gb = F[1:12, 9:16] = GridLayout()
    gc = F[1:12, 17:24] = GridLayout()
    
    gl1 = F[1,0] = GridLayout()
    gl2 = F[2,0] = GridLayout()
    gl3 = F[7,0] = GridLayout()

    gt1 = F[0,2] = GridLayout()
    gt2 = F[0,10] = GridLayout()
    gt3 = F[0,18] = GridLayout()

    gl = F[8:12,0] = GridLayout()

    #measurement
    #p = local_parameters(W,16,17,N,startyear)
    #mode_figure(ga,p,rich("TS [deg C]"))
    p = local_parameters(W,1,9,N,startyear)
    mode_figure(ga,p,rich("GPP [gC m", superscript("-2")," d",superscript("-1"),"]"))


    #resolved
    p = local_parameters(W,9,10,N,startyear)
    mode_figure(gb,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"))
    
    #unresolved
    p = local_parameters(W,7,3,N,startyear)
    mode_figure(gc,p,rich("RECO [gC m", superscript("-2")," d",superscript("-1"),"]"))

    bigfs = 26
    label1 = "signal and\n reconstruction"
    label2 = "relative power"
    label3 = "individual modes"

    title1 = "(a) measurement TS IT-Lav ENF"
    title2 = "(b) resolved NEE DE-Tha ENF"
    title3 = "(c) unresolved RECO CH-Dav ENF"

    for (label, layout) in zip([title1,title2,title3], [gt1,gt2,gt3])
        Label(layout[1, 1], label,
            fontsize = bigfs,
            font = :bold,
            #padding = (0, 5, 5, 0),
            #halign = :right
            )
    end

    for (label, layout) in zip([label1, label2, label3], [gl1, gl2, gl3])
        Label(layout[1, 1,], label,
            fontsize = bigfs,
            font = :bold,
            #padding = (0, 5, 5, 0),
            #halign = :right,
            #valign = :bottom,
            rotation=pi/2)
    end

    elem_1 = [LineElement(color = :grey68, linestyle = :solid)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid)]

    Legend(gl[1,1],
    [elem_1, elem_2, elem_3],
    ["Signal", "SSA", "NLSA"],
    labelsize = bigfs)

    for g in [ga, gb, gc]
        g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end

    for g in [gt1,gt2,gt3, gl1,gl2,gl3,gl]
        colgap!(g, 0)
        rowgap!(g, 0)
        g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end

    #colsize!(F.layout, 0, 15)
    #colsize!(F.layout, 2, 75)
    #colsize!(F.layout, 3, 75)
    #colsize!(F.layout, 4, 75)

    save(dir*savedirname*".png",F)
end

"""
outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass6/"
large_mode_figure("examples_lowpass6a")

outdir = "/net/scratch/lschulz/fluxfullset_midwithnee/"
large_mode_figure("examples_unfiltered")

"""



"""
combined heatmap
"""

#this doesnt work inside function scope, needs to be performed globally
function gather_trends()
    outdir = "/net/scratch/lschulz/fluxfullset_midwithnee/"
    tensors_raw = create_trend_tensors(N,W)

    outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass3/"
    tensors_3 = create_trend_tensors(N,W)

    outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass4/"
    tensors_4 = create_trend_tensors(N,W)

    outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass6/"
    tensors_6 = create_trend_tensors(N,W)

    jldsave(dir*"trend_tensors.jld2",
    tensors_raw = tensors_raw,
    tensors_3 = tensors_3,
    tensors_4 = tensors_4,
    tensors_6 = tensors_6
    )
end

function combined_heatmap()
    #ssa_h,nlsa_h,ssa_trends,nlsa_trends inside the tensors
    trends = load("/net/home/lschulz/logs/KW_2_20/trend_tensors.jld2")
    tensors_raw = trends["tensors_raw"]
    tensors_3 = trends["tensors_3"]
    tensors_4 = trends["tensors_4"]
    tensors_6 = trends["tensors_6"]


    F = Figure(resolution = (800,1200))
    
    g_a = F[1:3,1:10] = GridLayout()
    g_b = F[1:3,11:20] = GridLayout()
    g_c = F[4:6,1:10] = GridLayout()
    g_d = F[4:6,11:20] = GridLayout()
    g_e = F[7:9,1:10] = GridLayout()
    g_f = F[7:9,11:20] = GridLayout()
    g_g = F[10:12,1:10] = GridLayout()
    g_h = F[10:12,11:20] = GridLayout()

    g_cbar = F[0,1:20] = GridLayout()

    g_lab_1 = F[2,0] = GridLayout()
    g_lab_2 = F[5,0] = GridLayout()
    g_lab_3 = F[8,0] = GridLayout()
    g_lab_4 = F[11,0] = GridLayout()

    highclip = 8
    fs = 20
    fs_plot = 16


    function plot_heatmap(g,data)

        data = data[spots_list,vari_list]

        data[data .< 2] .= NaN

        x = 1:size(data,1)
        y = 1:size(data,2)

        ax,hm = heatmap(g[1,1],x,y,data,
            colormap = cgrad(:heat, highclip, categorical = true),
            colorrange = (1,highclip),
            highclip = :red,
            lowclip = :white,
            #colorrange = (minimum(valuetable),maximum(valuetable)),
            #aspect_ratio = 1,
            grid = true,
            framevisible = false,
            axis = (
                xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
                xticklabelrotation = pi/2,
                #xlabel="site",
                #ylabel = "variable",
                xlabelsize = fs,
                xticklabelsize = fs-4,
                yticks = (1:length(variables_names),variables_names),
                yticklabelsize = fs-4,
                ylabelsize = fs,
                #title = "SSA",
                #titlesize = fs,
            ))

        for i in 1:size(data)[1]
            str = data[i,:]
            map(enumerate(str)) do (j, c)
                if !isnan(c)
                    text!(ax,i,j,text="$(Int(c))",
                    align = (:center,:center),
                    color = :black,fontsize=fs_plot,
                    )
                end

            end
        end

        hidespines!(ax, :t,:r)
        return ax,hm

    end

    function select_variables()

        variables_names = ["GPP_DAY_1","GPP_DAY_2","GPP_NIGHT_1","GPP_NIGHT_2",
        "RECO_DAY_1","RECO_DAY_2","RECO_NIGHT_1","RECO_NIGHT_2",
        "NEE","NEE_NIGHT","NEE_DAY",
        "TA","SW_IN","LW_IN","SWC","TS"]

        #medium length
        spotslist = [
            "BE-Lon",
            "BE-Vie",
            "CH-Dav",
            "CH-Fru",
            "CH-Lae",
            "CH-Oe2",
            "DE-Geb",
            "DE-Gri",
            "DE-Hai",
            "DE-Tha",
            "DK-Sor",
            "ES-LJu",
            "FI-Hyy",
            "FR-Aur",
            "FR-Lam",
            "IT-BCi",
            "IT-Lav",
            "RU-Fyo",
        ]

        IGBP_list = [
            "CRO",
            "MF",
            "ENF",
            "GRA",
            "MF",
            "CRO",
            "CRO",
            "GRA",
            "DBF",
            "ENF",
            "DBF",
            "OSH",
            "ENF",
            "CRO",
            "CRO",
            "CRO",
            "ENF",
            "ENF",
        ]
        
        forest_list,grass_list = mask_IGBP(IGBP_list)
        spots_list = forest_list
        spotslist = spotslist[spots_list]
        IGBP_list = IGBP_list[forest_list]
        vari_list,variables_names = mask_vari(variables_names)

        return spots_list,spotslist,IGBP_list,vari_list,variables_names
    end

    spots_list,spotslist,IGBP_list,vari_list,variables_names = select_variables()

    ax_a,hm_a = plot_heatmap(g_a,tensors_raw[1])
    ax_b,hm_b = plot_heatmap(g_b,tensors_raw[2])
    ax_c,hm_c = plot_heatmap(g_c,tensors_6[1])
    ax_d,hm_d = plot_heatmap(g_d,tensors_6[2])
    ax_e,hm_e = plot_heatmap(g_e,tensors_4[1])
    ax_f,hm_f = plot_heatmap(g_f,tensors_4[2])
    ax_g,hm_g = plot_heatmap(g_g,tensors_3[1])
    ax_h,hm_h = plot_heatmap(g_h,tensors_3[2])

    #on the right, top3
    for ax in [ax_b,ax_d,ax_f]
        hidedecorations!(ax,grid = false,ticks=false)
    end

    #on the left, top 3
    for ax in [ax_a,ax_c,ax_e]
        hidexdecorations!(ax,grid = false,ticks=false)
    end

    #bottom figures
    hideydecorations!(ax_h,grid = false,ticks=false)

    for g in [g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_h,g_cbar,g_lab_1,g_lab_2,g_lab_3,g_lab_4]
        colgap!(g, 0)
        rowgap!(g, 0)
        #g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end


    Label(g_cbar[2,1],"SSA",
        fontsize = fs,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    Label(g_cbar[2,4],"NLSA",
        fontsize = fs,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)

    
    side_labels = ["unfiltered","lowpass [6/a]","lowpass [4/a]","lowpass [3/a]"]
    for (label,g) in zip(
        side_labels,[g_lab_1,g_lab_2,g_lab_3,g_lab_4]
        )
        Label(g[1, 1], label,
            fontsize = fs,
            font = :bold,
            rotation = pi/2,
            #padding = (0, 5, 5, 0),
            #halign = :right
            )
    end


    #Colorbar(g_cbar[1,1:6], hm_a, vertical = false,
    #label = "Harmonics",labelsize = fs
    #)

    subfigure_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
    gs = [g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_h]

    for (label, layout) in zip(subfigure_labels,gs)
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    for g in [g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_h]
        rowsize!(g, 1, Fixed(200))

    end

    save(dir*"heatmaps.png",F)
end

"""
seek regularity
"""
function regularity(signal)
    Ts = 1 / 365.25
    t0 = 0
    tmax = t0 + (N-1) * Ts
    t = t0:Ts:tmax
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    freq_domain_N = freqs[freqstart:freqend]
    years = ((1:N) ./ 365) .+ startyear

    spec_signal = abs.(fft(signal) |> fftshift)
    spec_signal /= norm(spec_signal)

    flist = [findall(x->(x-i)^2<0.001,freqs)[1] for i in 1:4]

    return spec_signal[flist],sum(spec_signal[freqend:end])
end

function heatmap_harmonics()
    spots = mask_IGBP(IGBP_list)[1]
    varis = mask_vari(variables_names)[1]
    harmonic_power = zeros(length(spots),length(varis),4)
    noise_power = zeros(length(spots),length(varis))
    for (i,spot) in enumerate(spots),(j,vari) in enumerate(varis)
        harmonic_power[i,j,:],noise_power[i,j] = regularity(wholedata[:,spot,vari])
    end

    F = Figure(resolution=(300,800))
    for i = 1:4
        ax = Axis(F[i,1])
        hm = heatmap!(ax,harmonic_power[:,:,i],
        colormap = cgrad(:heat, 100, categorical = true),
        #colorrange = (0,1),
        highclip = :red,
        lowclip = :white,
        #colorrange = (minimum(valuetable),maximum(valuetable)),
        #aspect_ratio = 1,
        grid = true,
        framevisible = false,
        )
        Colorbar(F[i,2], hm,
        vertical = true,
        #label = "Harmonics",labelsize = fs
        )
    end

    save(dir*"test.png",F)
end

function build_plots()


    spots = mask_IGBP(IGBP_list)[1]
    varis = mask_vari(variables_names)[1]
    spots = 1:18
    varis = 1:16
    power = zeros(length(spots),length(varis),5)


    for (i,spot) in enumerate(spots),(j,vari) in enumerate(varis)
        power[i,j,1:4],power[i,j,5] = regularity(wholedata[:,spot,vari])
    end

    trends = load(dir*"trend_tensors.jld2")
    tensors_raw = trends["tensors_raw"][1][spots,varis]
    tensors_3 = trends["tensors_3"][1][spots,varis]
    tensors_4 = trends["tensors_4"][1][spots,varis]
    tensors_6 = trends["tensors_6"][1][spots,varis]
    tensors = cat(tensors_raw,tensors_3,tensors_4,tensors_6,dims=3)

    markers = [
    :circle,
    :circle,
    :circle,
    :circle,
    :rect,
    :rect,
    :rect,
    :rect,
    :diamond,
    :diamond,
    :diamond,
    :utriangle,
    :dtriangle,
    :ltriangle,
    :rtriangle,
    :hexagon,

    ]

    F = Figure(resolution=(1200,1200))

    for (i,tit) in enumerate(["h1","h2","h3","h4","noise"]), (j,fil) in enumerate(["raw","low 3/a","low 4/a","low 6/a"])
        ax = Axis(F[i,j],title = fil,xlabel="intensity "*tit,ylabel="#H")
        scatter!(ax,power[:,:,i][:],tensors[:,:,j][:],
        marker = reshape(repeat(markers,18),16,18)[:],
        label = reshape(repeat(variables_names,18),16,18)[:])
    end

    save(dir*"SSA_H.png",F)


    #harmonic           NEES NEW POWER INDEX
    F = Figure()
    ax = Axis(F[1,1])
    hm = heatmap!(ax,
        sum(harmonic_power,dims=3)[:,:],
        colormap = cgrad(:heat, 100, categorical = true),
        #colorrange = (0,1),
        highclip = :red,
        lowclip = :white,
        #colorrange = (minimum(valuetable),maximum(valuetable)),
        #aspect_ratio = 1,
        grid = true,
        framevisible = false,
        )

    Colorbar(F[1,2], hm)

    #noise              NEES NEW POWER INDEX
    ax = Axis(F[2,1])
    hm = heatmap!(ax,
        noise_power,
        colormap = cgrad(:heat, 100, categorical = true),
        #colorrange = (0,1),
        highclip = :red,
        lowclip = :white,
        #colorrange = (minimum(valuetable),maximum(valuetable)),
        #aspect_ratio = 1,
        grid = true,
        framevisible = false,
        )

    Colorbar(F[2,2], hm)

    save(dir*"test.png",F)



end

"""
lag figure
needs to be run globally again (something about the p...?)
"""
function lag()
    varis = [16,1]
    var1 = "TS"
    var2 = "GPP"
    lags = -30:1:30

    for spot = mask_IGBP(IGBP_list)[1]
        try
            outdir = "/net/scratch/lschulz/fluxfullset_midwithnee/"
            storage_raw = zeros(N,length(varis),5)
            for (i,vari) in enumerate(varis)
                p = local_parameters(W,vari,spot,N,startyear)
                spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
                storage_raw[:,i,1] = signal
                storage_raw[:,i,2] = ssa_rec
                storage_raw[:,i,3] = nlsa_rec
                storage_raw[:,i,4] = ssa_trend_harm
                storage_raw[:,i,5] = nlsa_trend_harm
            end

            outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass6/"
            storage_f = zeros(N,length(varis),5)
            for (i,vari) in enumerate(varis)
                p = local_parameters(W,vari,spot,N,startyear)
                spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
                storage_f[:,i,1] = signal
                storage_f[:,i,2] = ssa_rec
                storage_f[:,i,3] = nlsa_rec
                storage_f[:,i,4] = ssa_trend_harm
                storage_f[:,i,5] = nlsa_trend_harm
            end

            savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_raw.jld2"
            raw_data = load(savedirname)["data"]
            savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_filtered.jld2"
            f_data = load(savedirname)["data_4a"]
            s_raw = cat([centralizer(raw_data[:,spot,j]) .+ i*2 for (i,j) in enumerate(varis)]...,dims=2)
            s_f = cat([centralizer(f_data[:,spot,j]) .+ i*2 for (i,j) in enumerate(varis)]...,dims=2)

            F = Figure(resolution=(1800,800))
            axis_l = [Axis(F[i,1:3]) for i =1:6]
            axis_r = [Axis(F[i,4]) for (j,i) = enumerate([1:3,4:6])]

            y = [s_raw[:,:]',
                storage_raw[:,:,4]',
                storage_raw[:,:,5]',
                s_f[:,:]',
                storage_f[:,:,4]',
                storage_f[:,:,5]']

            labels = ["signal raw","SSA trend raw","NLSA trend raw","signal f","SSA trend f","NLSA trend f"]

            for i = 1:6
                series!(axis_l[i],y[i],labels=[var1,var2])
                text!(axis_l[i],0.1,0.8,text=labels[i],space=:relative)
                hidedecorations!(axis_l[i],grid=false,ticks=false)
            end

            axislegend(axis_l[1],position=:rb)

            y_crosscor = [crosscor(y[i][1,:],y[i][2,:],lags) for i = 1:6]

            labels=["crosscorr raw","crosscorr lowpass 4/a"]

            for (i,j) = enumerate([1:3,4:6])
                m = cat(y_crosscor[j]...,dims=2)'
                series!(axis_r[i],lags,m,labels=["signal","SSA","NLSA"],color=[:black,:darkgreen,:purple])
                text!(axis_r[i],0.1,0.8,text=labels[i],space=:relative)
                maxi = [argmax(m,dims=2)[:,1][i][2] for i=1:3]
                scatter!(axis_r[i],Array(lags)[maxi],[m[i,maxi[i]] for i in 1:3],markersize=10)
                vlines!(axis_r[i],[Array(lags)[maxi[i]] for i in 1:3],ymax = 0.1,)
                axislegend(axis_r[i],position=:rb)
            end

            Label(F[0,1:4], "lag crosscorr for $(spotslist[spot]) $(IGBP_list[spot])
            $(variables_names[varis[1]]) and $(variables_names[varis[2]])",
                fontsize = 26,
                font = :bold,)
            save(dir*"cc/$(var1)_$(var2)/$(var1)_$(var2)_$(spot).png",F)
        catch
            println(spot)
        end
    end
end


function cc_maxima()
    varis = [16,9]
    lags = -60:1:60

    storage_lags = zeros(length(mask_IGBP(IGBP_list)[1]),6)
    storage_values = zeros(length(mask_IGBP(IGBP_list)[1]),6)

    for (i,spot) = enumerate(mask_IGBP(IGBP_list)[1])
        try
            outdir = "/net/scratch/lschulz/fluxfullset_midwithnee/"
            storage_raw = zeros(N,length(varis),5)
            for (i,vari) in enumerate(varis)
                p = local_parameters(W,vari,spot,N,startyear)
                spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
                storage_raw[:,i,1] = signal
                storage_raw[:,i,2] = ssa_rec
                storage_raw[:,i,3] = nlsa_rec
                storage_raw[:,i,4] = ssa_trend_harm
                storage_raw[:,i,5] = nlsa_trend_harm
            end

            outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass6/"
            storage_f = zeros(N,length(varis),5)
            for (i,vari) in enumerate(varis)
                p = local_parameters(W,vari,spot,N,startyear)
                spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
                storage_f[:,i,1] = signal
                storage_f[:,i,2] = ssa_rec
                storage_f[:,i,3] = nlsa_rec
                storage_f[:,i,4] = ssa_trend_harm
                storage_f[:,i,5] = nlsa_trend_harm
            end

            savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_raw.jld2"
            raw_data = load(savedirname)["data"]
            savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_filtered.jld2"
            f_data = load(savedirname)["data_4a"]
            s_raw = cat([centralizer(raw_data[:,spot,j]) .+ i*2 for (i,j) in enumerate(varis)]...,dims=2)
            s_f = cat([centralizer(f_data[:,spot,j]) .+ i*2 for (i,j) in enumerate(varis)]...,dims=2)

            y = [s_raw[:,:]',
            storage_raw[:,:,4]',
            storage_raw[:,:,5]',
            s_f[:,:]',
            storage_f[:,:,4]',
            storage_f[:,:,5]']

            y_crosscor = [crosscor(y[i][1,:],y[i][2,:],lags) for i = 1:6]


            for j = 1:6
                m = y_crosscor[j]
                maxi = argmax(m)
                storage_lags[i,j] = lags[maxi]
                storage_values[i,j] = m[maxi]
            end
        catch
        end
    end


    F = Figure(resolution=(800,800))
    axis = [Axis(F[i,1]) for i in 1:2]
    labels = ["raw","lowpass 4/a"]
    colors = [:black,:darkgreen,:purple]
    sites = spotslist[mask_IGBP(IGBP_list)[1]]
    for (i,j) = enumerate([1:3,4:6])
        points = Point2f.(storage_lags[:,j], storage_values[:,j])
        for k = 1:3
            scatter!(axis[i],points[:,k],markersize=10,
            color = colors[k],marker = :x)
            text!(axis[i],points[:,k],text=sites,
            color = colors[k],fontsize=12)
        end

        text!(axis[i],0.1,0.1,text=labels[i],space=:relative)
    end

    for (i,j) = enumerate([1:3,4:6])
        points = Point2f.(storage_lags[:,j], storage_values[:,j])
        for k = 1:3
            scatter!(axis[i],points[:,k],markersize=10,
            color = colors[k],marker = :x)
            text!(axis[i],points[:,k],text=sites,
            color = colors[k],fontsize=12)
        end

        text!(axis[i],0.1,0.1,text=labels[i],space=:relative)
    end


    save(dir*"TS_NEE.png",F)

end


"""
motivation plot paper
"""
function motivation_plot()
    #hainichen NEE
    p = local_parameters(W,9,9,N,startyear)

    ssa_trend_rc,nlsa_trend_rc = individual_harmonics(p,false)
    p = rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
    
    m = means[spot,vari]
    s = stds[spot,vari]

    ssa_fund = back_trafo(sum(ssa_trend_rc[:,1:2],dims=2)[:],m,s)
    nlsa_fund = back_trafo(sum(nlsa_trend_rc[:,1:2],dims=2)[:],m,s)

    seasonalitystring = ""
    seasonalitystring *= "ssa f \t"*string(round.(freq_ssa,digits=1))*"\n"
    seasonalitystring *= "nlsa f \t"*string(round.(freq_nlsa,digits=1))*"\n"
    seasonalitystring *= "ssa var \t"*string(round.(ssa_harm_var,digits=1))*"\n"
    seasonalitystring *= "nlsa var \t"*string(round.(nlsa_harm_var,digits=1))*"\n"

    #the figure

    offset = 15
    lw = 3
    lw_s = 1
    ms = 4
    fontsize = 22

    #custom yticks 2* -5:5 but with the absolute offset!

    tiks = Array(-2:2:2)
    yticks = (cat(tiks,tiks .+ offset,dims=1),
    cat(string.(tiks),string.(tiks),dims=1))
    F = Figure(resolution=(800,400))

    ax_time = Axis(F[1,1],
    xticks = Int.(floor.(years[1]:3:years[end])),
    yticks=yticks,
    yminorgridvisible = true,
    yminorgridstyle = :dot,
    yminorticks = IntervalsBetween(5),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel = "time (a)",
    ylabel = rich("NEE [gC m",superscript("-2"),"d",superscript("-1"),"]"),
    xlabelsize = fontsize,
    ylabelsize = fontsize,
    xticklabelsize = fontsize-4,
    yticklabelsize = fontsize-4,
    )


    """
    Harmonic
    """
    #signal background

    scatter!(ax_time,years,signal .+ offset,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms)
    li = lines!(ax_time,years,signal .+ offset,linewidth=lw_s,
    color=color_signal,label="signal")

    #ssa
    li = lines!(ax_time,years,ssa_trend_harm .+ offset,linewidth=lw,
    color=color_ssa,linestyle=:solid,label="SSA")

    #nlsa
    li = lines!(ax_time,years,nlsa_trend_harm .+ offset,linewidth=lw,
    color=color_nlsa,linestyle=:solid,label="NLSA")

    """
    Fund
    """
    scatter!(ax_time,years,signal,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms)
    lines!(ax_time,years,signal,linewidth=lw_s,
    color=color_signal)

    #ssa
    lines!(ax_time,years,ssa_fund,linewidth=lw,
    color=color_ssa,linestyle=:solid)


    #nlsa
    lines!(ax_time,years,nlsa_fund,linewidth=lw,
    color=color_nlsa,linestyle=:solid)


    text!(
        ax_time, 0, 1,
        text = "seasonal cycle", 
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = fontsize
    )
    text!(
        ax_time, 0, 0.02,
        text = "only fundamental", 
        align = (:left, :bottom),
        offset = (4, -2),
        space = :relative,
        fontsize = fontsize
    )

    text!(
        ax_time, 0.95, 1,
        text = "$(spotslist[spot]) $(igbpclass) $(varname)", 
        align = (:right, :top),
        #offset = (4, -2),
        space = :relative,
        fontsize = fontsize
    )

    #hideydecorations!(ax_time, label = false)
    hidespines!(ax_time, :t, :r)

    #tightlimits!(ax_time,Left(),Right())

    axislegend(ax_time,
        tellheight = false,
        tellwidth = false,
        margin = (4, 4, 4, 4),
        halign = :right, valign = :bottom, orientation = :horizontal,
        labelsize = fontsize-2,
    )
    save(dir*"motivation.png",F)

end


"""
QF "/net/scratch/lschulz/QF/qualityflags.jld2"
#flags = load("/net/scratch/lschulz/QF/qualityflags.jld2")["flags"]
"""

function flagstuff()
    """
    flags=wholedata,vars = QC_list
    flags
    "TA"
    "SW"
    "LW"
    "USTAR"
    "NEE"
    "GPP"
    "RECO"
    "SWC"
    "TS"
    varis
    "GPP_DAY_1"
    "GPP_DAY_2"
    "GPP_NIGHT_1"
    "GPP_NIGHT_2"
    "RECO_DAY_1"
    "RECO_DAY_2"
    "RECO_NIGHT_1"
    "RECO_NIGHT_2"
    5.57955  4.64126   0.179796    2.45708  4.53793     6.7842   8.41551  10.996   1.34751  1.71513

    julia> f3_harm_p,f3_noise_p = tensor_based_analytics(wholedata)
    "NEE"
    "NEE_NIGHT"
    "NEE_DAY"
    "TA"
    "SW_IN"
    "LW_IN"
    "SWC"
    "TS"

    """

    F = Figure()
    ax = Axis(F[1,1],limits=(1,N,0,25))
    spot = 3
    vari = 7
    vari_flag = 7
    lines!(ax,wholedata[:,spot,vari])
    for i = 1:9
        lines!(ax,flags[:,spot,i] .+ i)
    end
    save(dir*"flags.png",F)


    ZZ = zeros(18,6)
    skew = zeros(18,6)
    kurt = zeros(18,6)
    devi = zeros(18,6)
    for i = 1:18
        #"GPP"
        ZZ[i,1] = sum(flags[:,i,6])/N
        skew[i,1] = skewness(flags[:,i,6])
        kurt[i,1] = kurtosis(flags[:,i,6])
        devi[i,1] = std(flags[:,i,6])
        #"RECO"
        ZZ[i,2] = sum(flags[:,i,7])/N
        skew[i,2] = skewness(flags[:,i,7])
        kurt[i,2] = kurtosis(flags[:,i,7])
        devi[i,2] = std(flags[:,i,7])
        #"NEE"
        ZZ[i,3] = sum(flags[:,i,5])/N
        skew[i,3] = skewness(flags[:,i,5])
        kurt[i,3] = kurtosis(flags[:,i,5])
        devi[i,3] = std(flags[:,i,5])
        #"SW_IN"
        ZZ[i,4] = sum(flags[:,i,2])/N
        skew[i,4] = skewness(flags[:,i,2])
        kurt[i,4] = kurtosis(flags[:,i,2])
        devi[i,4] = std(flags[:,i,2])
        #"TS"
        ZZ[i,5] = sum(flags[:,i,9])/N
        skew[i,5] = skewness(flags[:,i,9])
        kurt[i,5] = kurtosis(flags[:,i,9])
        devi[i,5] = std(flags[:,i,9])
        #"SWC"
        ZZ[i,6] = sum(flags[:,i,8])/N
        skew[i,6] = skewness(flags[:,i,8])
        kurt[i,6] = kurtosis(flags[:,i,8])
        devi[i,6] = std(flags[:,i,8])
    end

    F = Figure(resolution=(600,800))
    for (i,sm) = enumerate([ZZ,devi,skew,kurt])
        ax = Axis(F[i,1])
        hm = heatmap!(ax,sm[mask_IGBP(IGBP_list)[1],:],
        colormap = cgrad(:heat, 100, categorical = true),
        #colorrange = (0,1),
        highclip = :red,
        lowclip = :white,
        #colorrange = (minimum(valuetable),maximum(valuetable)),
        #aspect_ratio = 1,
        grid = true,
        framevisible = false,
        )
        Colorbar(F[i,2], hm, label = ["average","deviation","skewness","kurtosis"][i])

    end
    save(dir*"flags.png",F)
end

"""
the tensors are in the dir * trend_tensors.jld2 file
run globally!
"""

function all_the tensorz()

    function long_deviation(vec)
        st = std(vec)
        me = mean(vec)
        window_length = 182
        Ar = embed_lag(Float32.(vec),window_length)
        for i in 1:size(Ar)[2]
            win = Ar[:,i]
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
        freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

        freqstart = findall(x->x>=1/12,freqs)[1]
        freqend = findall(x->x>=6,freqs)[1]

        spec_signal = abs.(fft(signal) |> fftshift)
        h_freqs = [findall(x->(x-i)^2<0.001,freqs)[1] for i in 1:4]
        harmonic_power = sum(spec_signal[h_freqs])
        noise_power = sum(spec_signal[freqend:end])
        power_sum = sum(spec_signal)

        return harmonic_power/power_sum,noise_power/power_sum
    end

    #signal tensor based analytics
    function tensor_based_analytics(tensor)
        N,n_spots,n_varis = size(tensor)

        tensors = [intra_regularity(tensor[:,i,j]) for i =1:n_spots,j=1:n_varis] 
        harmonic_power = [tensors[i,j][1] for i in 1:n_spots,j=1:n_varis]
        noise_power = [tensors[i,j][2] for i in 1:n_spots,j=1:n_varis]

        return harmonic_power,noise_power
    end

    #projection binning with few bins: for different heatmap symbols
    function projection_binning_intra(matrices)
        function doit(matrix)
            n_spots,n_varis = size(matrix)
            bins = 3
            maxi = maximum(matrix) *1.01
            mini = minimum(matrix)
            binsize = (maxi-mini)/bins
            bin_thresholds = [mini+i*binsize for i in 0:bins]

            binning = zeros(Int64,n_spots,n_varis)
            for i = 1:n_spots,j=1:n_varis
                w = fit(Histogram,[matrix[i,j],],bin_thresholds).weights
                #println(w)
                binning[i,j] = sum(w .* (1:bins))
            end
            return Int64.(binning)
        end
        return doit.(matrices)
    end

    #intra-comparison is not suitable for filtered noise
    #build comparison across all filters
    function projection_binning_inter(matrices)
        n_spots,n_varis = size(matrices[1])
        n_matrices=size(matrices)[1]
        bins = 3
        maxi = maximum(maximum.(matrices)) *1.01
        mini = minimum(minimum.(matrices))
        binsize = (maxi-mini)/bins
        bin_thresholds = [mini+i*binsize for i in 0:bins]

        binning = zeros(Int64,n_matrices,n_spots,n_varis)
        for m = 1:n_matrices,i = 1:n_spots,j=1:n_varis
            w = fit(Histogram,[matrices[m][i,j],],bin_thresholds).weights
            #println(w)
            binning[m,i,j] = sum(w .* (1:bins))
        end
        return [binning[m,:,:] for m = 1:n_matrices]
    end


    savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/"*"fluxdata_raw.jld2"
    wholedata = SharedArray{Float32}(load(savedirname)["data"])
    raw_harm_p,raw_noise_p = tensor_based_analytics(wholedata)

    savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_filtered.jld2"
    wholedata_f4 = SharedArray{Float32}(load(savedirname)["data_4a"])
    f4_harm_p,f4_noise_p = tensor_based_analytics(wholedata)
    wholedata_f3 = SharedArray{Float32}(load(savedirname)["data_3a"])
    f3_harm_p,f3_noise_p = tensor_based_analytics(wholedata)
    wholedata_f6 = SharedArray{Float32}(load(savedirname)["data_6a"])
    f6_harm_p,f6_noise_p = tensor_based_analytics(wholedata)

    flags = load("/net/scratch/lschulz/QF/qualityflags.jld2")["flags"]
    flags[flags .< 0] .= 0
    flags[flags .> 1] .= 1
    lo = [long_deviation(flags[:,i,j]) for i = 1:18,j = 1:9]
    #lo[spots,[6,7,5,2,9,8]]

    noises = [raw_noise_p[spots,vars],f3_noise_p[spots,vars],f4_noise_p[spots,vars],f6_noise_p[spots,vars]]
    harms = [raw_harm_p[spots,vars],f3_harm_p[spots,vars],f4_harm_p[spots,vars],f6_harm_p[spots,vars]]
    noises_scaled = projection_binning_inter(noises)
    harms_scaled = projection_binning_inter(harms)

    jldsave(dir*"intra_tensor_based_analytics.jld2",
        raw_harm_p = harms_scaled[1],
        f3_harm_p = harms_scaled[2],
        f4_harm_p = harms_scaled[3],
        f6_harm_p = harms_scaled[4],
        raw_noise_p = noises_scaled[1],
        f3_noise_p = noises_scaled[2],
        f4_noise_p = noises_scaled[3],
        f6_noise_p = noises_scaled[4],
        lo = lo[spots,[6,7,5,2,9,8]],
    )


    raw_entropy = [sampen3(data_raw[:,i,j]) for i = 1:18,j = 1:16]
    f3_entropy = [sampen3(data_f3[:,i,j]) for i = 1:18,j = 1:16]
    f4_entropy = [sampen3(data_f4[:,i,j]) for i = 1:18,j = 1:16]
    f6_entropy = [sampen3(data_f6[:,i,j]) for i = 1:18,j = 1:16]
    
    
    
    data = load("/net/scratch/lschulz/data/time_series.jld2")
    data_raw = data["data_raw"]
    data_f3 = data["data_f3" ]
    data_f4 = data["data_f4"]
    data_f6 = data["data_f6"]

    raw_harm_p,raw_noise_p = tensor_based_analytics(data_raw)

    f4_harm_p,f4_noise_p = tensor_based_analytics(data_f4)

    f3_harm_p,f3_noise_p = tensor_based_analytics(data_f3)

    f6_harm_p,f6_noise_p = tensor_based_analytics(data_f6)

    noises = [raw_noise_p[spots,vars],f3_noise_p[spots,vars],f4_noise_p[spots,vars],f6_noise_p[spots,vars]]
    harms = [raw_harm_p[spots,vars],f3_harm_p[spots,vars],f4_harm_p[spots,vars],f6_harm_p[spots,vars]]
    entropies = [raw_entropy[spots,vars],f3_entropy[spots,vars],f4_entropy[spots,vars],f6_entropy[spots,vars]]
    noises_scaled = projection_binning_inter(noises)
    harms_scaled = projection_binning_inter(harms)
    entropies_scaled = projection_binning_inter(entropies)

    jldsave("/net/scratch/lschulz/data/data_characteristics.jld2",
        f3_harm_p = f3_harm_p,
        f4_harm_p = f4_harm_p,
        f6_harm_p = f6_harm_p,
        raw_harm_p = raw_harm_p,
        f3_noise_p = f3_noise_p,
        f4_noise_p = f4_noise_p,
        f6_noise_p = f6_noise_p,
        raw_noise_p = raw_noise_p,
        raw_entropy = raw_entropy,
        f3_entropy = f3_entropy,
        f4_entropy = f4_entropy,
        f6_entropy = f6_entropy,
        raw_harm_p_scaled = harms_scaled[1],
        f3_harm_p_scaled = harms_scaled[2],
        f4_harm_p_scaled = harms_scaled[3],
        f6_harm_p_scaled = harms_scaled[4],
        raw_noise_p_scaled = noises_scaled[1],
        f3_noise_p_scaled = noises_scaled[2],
        f4_noise_p_scaled = noises_scaled[3],
        f6_noise_p_scaled = noises_scaled[4],
        raw_entropy_scaled = entropies_scaled[1],
        f3_entropy_scaled = entropies_scaled[2],
        f4_entropy_scaled = entropies_scaled[3],
        f6_entropy_scaled = entropies_scaled[4],
        artifacts = artifacts,
    )

    function rle_encode(match_matrix)
        N = size(match_matrix, 1)
        encoded = Vector{Int}(undef, N)
        for i in 1:N
            encoded[i] = sum(match_matrix[i, i+1:end])
        end
        return encoded
    end

    function calculate_match_matrix(L, m, r)
        N = length(L)
        match_matrix = falses(N - m + 1, N - m + 1)
        for i in 1:N - m
            for j in i + 1:N - m + 1
                if abs(L[i] - L[j]) ≤ r
                    match_matrix[i, j] = true
                else
                    break
                end
            end
        end
        return match_matrix
    end
    
    function rle_encode(match_matrix)
        N = size(match_matrix, 1)
        encoded = Vector{Int}(undef, N)
        for i in 1:N
            encoded[i] = sum(match_matrix[i, i + 1:N + 1])
        end
        return encoded
    end
    
    function sampen3(L)
        m = 2
        r = 2 * std(L)
        N = length(L)
        B = 0
        A = 0
    
        # Calculate distance matrix for B
        match_matrix_B = calculate_match_matrix(L, m, r)
        encoded_B = rle_encode(match_matrix_B)
        B = sum(encoded_B)
    
        # Calculate distance matrix for A
        m += 1
        match_matrix_A = calculate_match_matrix(L, m, r)
        encoded_A = rle_encode(match_matrix_A)
        A = sum(encoded_A)
    
        # Handle division by zero
        if A == 0 || B == 0
            return NaN
        end
    
        # Return SampEn
        return -log(A / B)
    end

    
end

"""
heatmaps with symbols
this is not working out at the moment
"""

function specified_combined_heatmap()
    #ssa_h,nlsa_h,ssa_trends,nlsa_trends inside the tensors
    trends = load("/net/home/lschulz/logs/KW_2_20/trend_tensors.jld2")
    tensors_raw = trends["tensors_raw"]
    tensors_3 = trends["tensors_3"]
    tensors_4 = trends["tensors_4"]
    tensors_6 = trends["tensors_6"]


    F = Figure(resolution = (800,1200))
    
    g_a = F[1:3,1:10] = GridLayout()
    g_b = F[1:3,11:20] = GridLayout()
    g_c = F[4:6,1:10] = GridLayout()
    g_d = F[4:6,11:20] = GridLayout()
    g_e = F[7:9,1:10] = GridLayout()
    g_f = F[7:9,11:20] = GridLayout()
    g_g = F[10:12,1:10] = GridLayout()
    g_h = F[10:12,11:20] = GridLayout()

    g_cbar = F[0,1:20] = GridLayout()

    g_lab_1 = F[2,0] = GridLayout()
    g_lab_2 = F[5,0] = GridLayout()
    g_lab_3 = F[8,0] = GridLayout()
    g_lab_4 = F[11,0] = GridLayout()

    highclip = 8
    fs = 20
    fs_plot = 16


    function plot_heatmap(g,data)

        data = data[spots_list,vari_list]

        data[data .< 2] .= NaN

        x = 1:size(data,1)
        y = 1:size(data,2)

        ax,hm = heatmap(g[1,1],x,y,data,
            colormap = cgrad(:heat, highclip, categorical = true),
            colorrange = (1,highclip),
            highclip = :red,
            lowclip = :white,
            #colorrange = (minimum(valuetable),maximum(valuetable)),
            #aspect_ratio = 1,
            grid = true,
            framevisible = false,
            axis = (
                xticks = (1:length(spotslist),spotslist.*" ".*IGBP_list),
                xticklabelrotation = pi/2,
                #xlabel="site",
                #ylabel = "variable",
                xlabelsize = fs,
                xticklabelsize = fs-4,
                yticks = (1:length(variables_names),variables_names),
                yticklabelsize = fs-4,
                ylabelsize = fs,
                #title = "SSA",
                #titlesize = fs,
            ))

        for i in 1:size(data)[1]
            str = data[i,:]
            map(enumerate(str)) do (j, c)
                if !isnan(c)
                    text!(ax,i,j,text="$(Int(c))",
                    align = (:center,:center),
                    color = :black,fontsize=fs_plot,
                    )
                end

            end
        end

        hidespines!(ax, :t,:r)
        return ax,hm

    end

    function select_variables()

        variables_names = ["GPP_DAY_1","GPP_DAY_2","GPP_NIGHT_1","GPP_NIGHT_2",
        "RECO_DAY_1","RECO_DAY_2","RECO_NIGHT_1","RECO_NIGHT_2",
        "NEE","NEE_NIGHT","NEE_DAY",
        "TA","SW_IN","LW_IN","SWC","TS"]

        #medium length
        spotslist = [
            "BE-Lon",
            "BE-Vie",
            "CH-Dav",
            "CH-Fru",
            "CH-Lae",
            "CH-Oe2",
            "DE-Geb",
            "DE-Gri",
            "DE-Hai",
            "DE-Tha",
            "DK-Sor",
            "ES-LJu",
            "FI-Hyy",
            "FR-Aur",
            "FR-Lam",
            "IT-BCi",
            "IT-Lav",
            "RU-Fyo",
        ]

        IGBP_list = [
            "CRO",
            "MF",
            "ENF",
            "GRA",
            "MF",
            "CRO",
            "CRO",
            "GRA",
            "DBF",
            "ENF",
            "DBF",
            "OSH",
            "ENF",
            "CRO",
            "CRO",
            "CRO",
            "ENF",
            "ENF",
        ]
        
        forest_list,grass_list = mask_IGBP(IGBP_list)
        spots_list = forest_list
        spotslist = spotslist[spots_list]
        IGBP_list = IGBP_list[forest_list]
        vari_list,variables_names = mask_vari(variables_names)

        return spots_list,spotslist,IGBP_list,vari_list,variables_names
    end

    function write_symbol(ax,data,symbol,x_shift,y_shift) #5 bins
        symbol_fontsize_list = 0:2:8
        for i in 1:size(data)[1]
            str = data[i,:]
            map(enumerate(str)) do (j, c)
            fontsize = symbol_fontsize_list[Int(c)]
            text!(ax,i+x_shift,j+y_shift,text=symbol,
            align = (:center,:center),
            color = :black,
            fontsize=fontsize,
            )

            end
        end
    end

    spots_list,spotslist,IGBP_list,vari_list,variables_names = select_variables()

    ax_a,hm_a = plot_heatmap(g_a,tensors_raw[1])
    ax_b,hm_b = plot_heatmap(g_b,tensors_raw[2])
    ax_c,hm_c = plot_heatmap(g_c,tensors_6[1])
    ax_d,hm_d = plot_heatmap(g_d,tensors_6[2])
    ax_e,hm_e = plot_heatmap(g_e,tensors_4[1])
    ax_f,hm_f = plot_heatmap(g_f,tensors_4[2])
    ax_g,hm_g = plot_heatmap(g_g,tensors_3[1])
    ax_h,hm_h = plot_heatmap(g_h,tensors_3[2])

    #on the right, top3
    for ax in [ax_b,ax_d,ax_f]
        hidedecorations!(ax,grid = false,ticks=false)
    end

    #on the left, top 3
    for ax in [ax_a,ax_c,ax_e]
        hidexdecorations!(ax,grid = false,ticks=false)
    end

    #bottom figures
    hideydecorations!(ax_h,grid = false,ticks=false)

    for g in [g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_h,g_cbar,g_lab_1,g_lab_2,g_lab_3,g_lab_4]
        colgap!(g, 0)
        rowgap!(g, 0)
        #g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end


    Label(g_cbar[2,1],"SSA",
        fontsize = fs,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    Label(g_cbar[2,4],"NLSA",
        fontsize = fs,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)

    
    side_labels = ["unfiltered","lowpass [6/a]","lowpass [4/a]","lowpass [3/a]"]
    for (label,g) in zip(
        side_labels,[g_lab_1,g_lab_2,g_lab_3,g_lab_4]
        )
        Label(g[1, 1], label,
            fontsize = fs,
            font = :bold,
            rotation = pi/2,
            #padding = (0, 5, 5, 0),
            #halign = :right
            )
    end


    #Colorbar(g_cbar[1,1:6], hm_a, vertical = false,
    #label = "Harmonics",labelsize = fs
    #)

    subfigure_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
    gs = [g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_h]

    for (label, layout) in zip(subfigure_labels,gs)
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    for g in [g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_h]
        rowsize!(g, 1, Fixed(200))

    end



    specifications = load(dir*"tensor_based_analytics.jld2")
    raw_harm_p = specifications["raw_harm_p"]
    raw_noise_p = specifications["raw_noise_p"]
    f3_harm_p = specifications["f3_harm_p"]
    f3_noise_p = specifications["f3_noise_p"]
    f4_harm_p = specifications["f4_harm_p"]
    f4_noise_p = specifications["f4_noise_p"]
    f6_harm_p = specifications["f6_harm_p"]
    f6_noise_p = specifications["f6_noise_p"]
    lo = specifications["lo"]

    filtered_specifications = [[raw_harm_p,raw_noise_p],
                        [f6_harm_p,f6_noise_p],
                        [f4_harm_p,f4_noise_p],
                        [f3_harm_p,f3_noise_p]]


    filtered_axis = [[ax_a,ax_b],
                    [ax_c,ax_d],
                    [ax_e,ax_f],
                    [ax_g,ax_h]]

    syms = [L"f",L"$\eta$"]
    x_shifts = [0.3,-0.3]
    y_shifts = [0.3,-0.3]

    
    for (axis,datas) = zip(filtered_axis,filtered_specifications)
        for (a,ax) = enumerate(axis), (d,data) = enumerate(datas)
            println("axis $a data $d value $(sum(datas))")
            symbol = syms[d]
            x_shift = x_shifts[d]
            y_shift = y_shifts[d]
            write_symbol(ax,data,symbol,x_shift,y_shift)
        end
    end

    save(dir*"heatmaps.png",F)
end



function characterization_map()

    function write_symbol(ax,data,symbol,x_shift,y_shift) #5 bins
        symbol_fontsize_list = [0,8,12,16,20]
        for i in 1:size(data)[1]
            str = data[i,:]
            map(enumerate(str)) do (j, c)
            fontsize = symbol_fontsize_list[Int(c)]
            text!(ax,i+x_shift,j+y_shift,text=symbol,
            align = (:center,:center),
            color = :black,
            fontsize=fontsize,
            )

            end
        end
    end

    function artifact_symbol(ax,data,x_shift,y_shift)
        symbol = L"X"
        fs = 20
        for i in 1:size(data)[1]
            str = data[i,:]
            map(enumerate(str)) do (j, c)
                if c
                    text!(ax,i+x_shift,j+y_shift,text=symbol,
                    align = (:center,:center),
                    color = :black,
                    fontsize=fs,
                    )
                end
            end
        end
    end

    specifications = load(dir*"intra_tensor_based_analytics.jld2")
    raw_harm_p = specifications["raw_harm_p"]
    raw_noise_p = specifications["raw_noise_p"]
    lo = specifications["lo"]
    samp_en = load(dir*"sample_entropy.jld2")["samp_en"]

    trends = load("/net/home/lschulz/logs/KW_2_20/trend_tensors.jld2")
    tensors_raw = trends["tensors_raw"]


    F = Figure()
    ax = Axis(F[1,1],xticks = (1:9,spotslist.*" ".*IGBP_list),
    xticklabelrotation = pi/2,
    yticks = (1:6,variables_names),
    yticklabelrotation = pi/2,
    )




    h = heatmap!(ax,tensors_raw[2][spots,vars],
        colormap = cgrad(:heat, highclip, categorical = true),
        colorrange = (1,highclip),
        highclip = :red,
        lowclip = :white,
        #colorrange = (minimum(valuetable),maximum(valuetable)),
        #aspect_ratio = 1,
        grid = true,
        framevisible = false,
    )

    Colorbar(F[1,2], h)
    write_symbol(ax,raw_harm_p,L"$\mathbb{H}$",0.3,0.3)
    write_symbol(ax,raw_noise_p,L"$\xi$",-0.3,-0.3)
    write_symbol(ax,samp_en,L"$E$",-0.3,0.3)
    artifact_symbol(ax,lo,0.3,-0.3)

    vlines!(ax,(0:9) .+ 0.5)
    hlines!(ax,(0:6) .+ 0.5)
    save(dir*"characterization_map.png",F)
end

function sampen3(L)
    m = 2
    r = 2 * std(L)
    N = length(L)
    B = 0
    A = 0

    # Calculate distance matrix for B
    diff_matrix_B = []
    for i in 1:N - m
        xmi = L[i : i + m - 1]
        diffs = []
        for j in i + 1:N - m + 1
            xmj = L[j : j + m - 1]
            diffs = [maximum(abs.(xmi .- xmj))]
        end
        diff_matrix_B = [diff_matrix_B; diffs]
    end
    match_matrix_B = [all(diff .≤ r) for diff in diff_matrix_B]
    B = sum(match_matrix_B)

    # Calculate distance matrix for A
    m += 1
    diff_matrix_A = []
    for i in 1:N - m
        xm = L[i : i + m - 1]
        diffs = []
        for j in i + 1:N - m + 1
            xj = L[j : j + m - 1]
            diffs = [maximum(abs.(xm .- xj))]
        end
        diff_matrix_A = [diff_matrix_A; diffs]
    end
    match_matrix_A = [all(diff .≤ r) for diff in diff_matrix_A]
    A = sum(match_matrix_A)

    # Handle division by zero
    if A == 0 || B == 0
        return NaN
    end

    # Return SampEn
    return -log(A / B)
end



"""
large figure with quality flags

outdir = "/net/scratch/lschulz/fluxfullset_midwithnee_lowpass6/"
large_mode_figure_flags("examples_lowpass6a")

outdir = "/net/scratch/lschulz/fluxfullset_midwithnee/"
large_mode_figure_flags("examples_unfiltered")

"""

function mode_figure_flags(F,p,varname_resolved,flags)

    p = rescale_local_parameters(p)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    modenumber = 12
    modenumber_k = 48
    xstep_k = 16
    k_xticks = (xstep_k:xstep_k:modenumber_k,string.(xstep_k:xstep_k:modenumber_k))
    spec_xticks = (1:6,string.(1:6))
    smallfs = 16
    fs = 23
    speclimits = (freq_domain_N[1],6.05,10^-4,1)
    klimits = (0,modenumber_k+1,10^-4,1)
    modeslimits = (0,2556/365.25,-0.04,(2*modenumber+2.5)*0.08)
    modeyticks = (vcat(0.08:0.08*2:(modenumber)*0.08,(modenumber+2)*0.08:0.08*2:(2*modenumber+1)*0.08),string.(vcat(1:2:modenumber,1:2:modenumber)))
    freqlimits = (1/10, 7,-0.2,modenumber*2+2.5)
    gapsize = 8
    lw = 3
    ms = 12
    scatterstep = 10
    color_signal = "grey68"

    

    #figure

    ax_time = Axis(F[1:2,1:4],
    xticks=Int.(floor.(years[1]:3:years[end])),
    limits = (years[1],years[end],minimum(signal)*1.1,maximum(signal)*1.1),
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time [a]",
    ylabel=varname_resolved,
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_time,ticks = false,ticklabels=false,grid=false,label=false)

    ax_spec = Axis(F[3:4,1:3],yscale=log10,
    limits = speclimits,
    xticks = spec_xticks,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(4),
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_spec,ticks = false,ticklabels=false,grid=false)

    ax_k = Axis(F[3:4,4],yscale=log10,
    limits = klimits,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(2),
    xlabel="dim",
    xlabelsize = fs,
    xticks = k_xticks,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_k,grid=false)


    ax_modes = Axis(F[5:9,1:2],#yticksvisible = false,
    #yticklabelsvisible = false,
    xticks=0:Int(floor(W/365.25)),
    yticks = modeyticks,
    limits = modeslimits,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(2),
    xlabel="time [a]",
    ylabel="individual modes",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_modes,ticks = false,ticklabels=false,grid=false)

    ax_freq = Axis(F[5:9,3:4],
    limits=freqlimits,#,yscale=log10
    xticks=1:7,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs,
    yticklabelsize = fs,
    ylabelsize = fs,)

    hideydecorations!(ax_freq)

    #ax_labels = Axis(F[3:4,5])
    #hidedecorations!(ax_labels)

    #flags
    vlines!(ax_time,years,color=flags,ymin = 0.0, ymax = 1.0,colormap = :grays)
    
    #plotting
    #scatter!(ax_time,years[1:scatterstep:end],signal[1:scatterstep:end],color=color_signal,markersize = ms,marker=:x)
    signal_l = lines!(ax_time,years,signal,color=color_signal,linewidth=lw,label = "signal")
    ssa_l = lines!(ax_time,years,ssa_trend_harm,color=color_ssa,linewidth=lw,label="SSA")
    nlsa_l = lines!(ax_time,years,nlsa_trend_harm,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain_N,spec_signal,color=color_signal,linewidth=lw,label = "signal")
    signal_s = scatter!(ax_spec,freq_domain_N,spec_signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_spec,freq_domain_N,spec_ssa_rc,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain_N,spec_nlsa_rc,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_k,1:modenumber_k,ssa_cap_var[1:modenumber_k],color=color_ssa,linewidth=lw)
    lines!(ax_k,1:modenumber_k,nlsa_cap_var[1:modenumber_k],color=color_nlsa,linewidth=lw)




    #modes 

    color_fit = "grey"
    for i = 1:modenumber

        mode = ssa_Eof[:,i]
        Four = spec_ssa_eof[:,i]
        #Four ./= maximum(abs.(Four)) 

        years_modes = (1:length(mode)) ./ 365.25


        lines!(ax_modes,years_modes,mode .+ (i*0.08),
        color=color_ssa)
        #lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_ssa[i]) .+ i,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i, #[freqstart:end]
        color=color_ssa)

        mode = nlsa_Eof[:,i] 
        Four = spec_nlsa_eof[:,i]


        lines!(ax_modes,years_modes,mode .+ ((modenumber +1+i)*0.08),
        color=color_nlsa)
        #lines!(ax_freq,freqs_w,gauss(freqs_w,gaussian_nlsa[i]) .+ i .+ modenumber .+1,color=color_fit,linewidth=5)
        lines!(ax_freq,freqs_w,abs.(Four) .+ i .+ modenumber .+1, #[freqstart:end]
        color=color_nlsa)

    end

    #height and freq
    boxstring = "f \t height \n"#"peaks \n\n f \t height \n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_nlsa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_nlsa[k][2],digits=1))
        boxstring *= "\n"
    end
    boxstring*= "\n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_ssa[k][1],digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_ssa[k][2],digits=1))
        boxstring *= "\n"
    end
    linkyaxes!(ax_spec, ax_k)

    #colgap!(F.layout, 1, gapsize) #colgap!(fig.layout, 1, Relative(0.15))
    #colgap!(F.layout, 2, gapsize)
    #colgap!(F.layout, 3, gapsize)
    #rowgap!(F.layout, 1, gapsize)
    #rowgap!(F.layout, 2, gapsize)
    #rowgap!(F.layout, 3, gapsize)

    return F
end

function large_mode_figure_flags(savedirname)
    F = Figure(resolution=(2400,1000))

    ga = F[1:12, 1:8] = GridLayout()
    gb = F[1:12, 9:16] = GridLayout()
    gc = F[1:12, 17:24] = GridLayout()
    
    gl1 = F[1,0] = GridLayout()
    gl2 = F[2,0] = GridLayout()
    gl3 = F[7,0] = GridLayout()

    gt1 = F[0,2] = GridLayout()
    gt2 = F[0,10] = GridLayout()
    gt3 = F[0,18] = GridLayout()

    gl = F[8:12,0] = GridLayout()

    flags = load("/net/scratch/lschulz/QF/qualityflags.jld2")["flags"]
    flags[flags .< 0] .= 0
    flags[flags .> 1] .= 1
    """
    local, filtered coordinates!
    """
    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]
    #measurement
    spot = 9
    var = 3
    p = local_parameters(W,vars[var],spots[spot],N,startyear)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(ga,p,rich("GPP [gC m", superscript("-2")," d",superscript("-1"),"]"),flag)
    title1 = "(a) $(variables_names[vars[var]]) $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])"

    #resolved
    spot = 2
    var = 3
    p = local_parameters(W,vars[var],spots[spot],N,startyear)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(gb,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag)
    title2 = "(b) $(variables_names[vars[var]]) $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])"

    #unresolved
    spot = 1
    var = 2
    p = local_parameters(W,vars[var],spots[spot],N,startyear)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(gc,p,rich("RECO [gC m", superscript("-2")," d",superscript("-1"),"]"),flag)
    title3 = "(c) $(variables_names[vars[var]]) $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])"

    bigfs = 26
    label1 = "signal and\n reconstruction"
    label2 = "relative power"
    label3 = "individual modes"



    for (label, layout) in zip([title1,title2,title3], [gt1,gt2,gt3])
        Label(layout[1, 1], label,
            fontsize = bigfs,
            font = :bold,
            #padding = (0, 5, 5, 0),
            #halign = :right
            )
    end

    for (label, layout) in zip([label1, label2, label3], [gl1, gl2, gl3])
        Label(layout[1, 1,], label,
            fontsize = bigfs,
            font = :bold,
            #padding = (0, 5, 5, 0),
            #halign = :right,
            #valign = :bottom,
            rotation=pi/2)
    end

    elem_1 = [LineElement(color = :grey68, linestyle = :solid)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid)]

    Legend(gl[1,1],
    [elem_1, elem_2, elem_3],
    ["Signal", "SSA", "NLSA"],
    labelsize = bigfs)

    for g in [ga, gb, gc]
        g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end

    for g in [gt1,gt2,gt3, gl1,gl2,gl3,gl]
        colgap!(g, 0)
        rowgap!(g, 0)
        g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end

    #colsize!(F.layout, 0, 15)
    #colsize!(F.layout, 2, 75)
    #colsize!(F.layout, 3, 75)
    #colsize!(F.layout, 4, 75)

    save(dir*savedirname*".png",F)
end

for spot=1:9,var = 1:6
    gb = Figure()
    p = local_parameters(W,vars[var],spots[spot],N,startyear)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(gb,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag)
    save(dir*"mode_figure_flags_filtered/$(spot)_$var.png",gb)
end


