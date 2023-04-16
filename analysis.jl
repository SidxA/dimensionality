include("/net/home/lschulz/dimensionality/fundamentals.jl")
using CairoMakie

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
    "unknown",
    "unknown",
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

outdir="/net/scratch/lschulz/fluxfullset_midwithnee/"
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

    return fit_gauss(freqs_w,spec)
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

normalizer(x) = x./maximum(x)

function local_parameters(W,vari,spot,N,startyear)

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
    ssa_cap_var = sum(ssa_lambda)
    ssa_rec = sum(ssa_RC,dims=2)[:]

    nlsa_lambda = file_nlsa["lambda"]
    nlsa_indices = sortperm(nlsa_lambda,rev=true)
    nlsa_Eof = file_nlsa["EOF"][:,nlsa_indices]
    nlsa_PC = file_nlsa["PC"][:,nlsa_indices]
    nlsa_RC = file_nlsa["RC"][:,nlsa_indices]
    nlsa_lambda = nlsa_lambda[nlsa_indices]
    nlsa_cap_var = sum(nlsa_lambda)
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
summarizing quantification
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


function heatmap_table(f,valuetable,titlestring)
    stepnr = 10
    stepsize = (maximum(valuetable)-minimum(valuetable))/(stepnr) 
    bins = LinRange(minimum(valuetable)-stepsize,maximum(valuetable)+stepsize,stepnr+2)

    cgradient = cgrad(:thermal, length(bins), categorical = true)
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

ssa_harmonics,nlsa_harmonics,ssa_max_slope,nlsa_max_slope,ssa_argmax_slope,nlsa_argmax_slope = create_valuetables(N,W)

dg = 3
valuetable_list_ssa = [
    Int.(ssa_harmonics),
    round.(ssa_max_slope,digits=dg),
    round.(ssa_argmax_slope,digits=dg),
]
valuetable_list_nlsa = [
    Int.(nlsa_harmonics),
    round.(nlsa_max_slope,digits=dg),
    round.(nlsa_argmax_slope,digits=dg),
]

value_title_list_ssa = [
    "SSA\n #harmonics",
    "SSA\nlin max slope",
    "SSA\nlin argmax slope",

]

value_title_list_nlsa = [
    "NLSA\n #harmonics",
    "NLSA\nlin max slope",
    "NLSA\nlin argmax slope",

]
fs = 12
padnr = 5

f = Figure(resolution=(800,1400))
for i = 1:3
    valuetable = valuetable_list_ssa[i]
    titlestring = value_title_list_ssa[i]
    heatmap_table(f[i,1],valuetable,titlestring)
end
save(dir*"ssa.png",f)
f = Figure(resolution=(800,1400))
for i = 1:3
    valuetable = valuetable_list_nlsa[i]
    titlestring = value_title_list_nlsa[i]
    heatmap_table(f[i,1],valuetable,titlestring)
end
save(dir*"nlsa.png",f)
