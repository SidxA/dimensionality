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

#the motivation figure
function motivation_plot(savedirname,outdir)
    #hainichen NEE
    spot = 9
    vari = 9
    p = local_parameters(spot,vari,data_tensor,outdir)

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

    back_trafo(data,mean,std) = (data.*std).+mean


    ssa_trend_rc,nlsa_trend_rc = individual_harmonics(p,false)
    p = rescale_local_parameters(p,data_tensor)

    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p
    
    m = mean(data_raw[:,spot,vari])
    s = std(data_raw[:,spot,vari])

    ssa_fund = back_trafo(sum(ssa_trend_rc[:,1:2],dims=2)[:],m,s)
    nlsa_fund = back_trafo(sum(nlsa_trend_rc[:,1:2],dims=2)[:],m,s)

    seasonalitystring = ""
    seasonalitystring *= "ssa f \t"*string(round.(freq_ssa,digits=1))*"\n"
    seasonalitystring *= "nlsa f \t"*string(round.(freq_nlsa,digits=1))*"\n"
    seasonalitystring *= "ssa var \t"*string(round.(ssa_harm_var,digits=1))*"\n"
    seasonalitystring *= "nlsa var \t"*string(round.(nlsa_harm_var,digits=1))*"\n"

    #the figure

    offset = 15

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
    save(savedirname,F)

end

"""
the following figures are built using the reduced number of spots and variables
"""

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

#individual overview plot showing the time series, the spectrum, the modes and the spectra of the modes
#need to be given a Figure() as F in order to combine them to a single figure
#varname resolved is the name of the variable together with the corresponding unit, making up the time series axis label
function mode_figure_flags(F,p,varname_resolved,flags,data_tensor)

    p = rescale_local_parameters(p,data_tensor)
    spot,W,vari,years,varname,igbpclass,freq_domain_N,freq_domain_w,freqs_w,freqs,signal,ssa_Eof,nlsa_Eof,nlsa_eps,ssa_rec,nlsa_rec,ssa_cap_var,nlsa_cap_var,spec_signal,spec_ssa_rc,spec_nlsa_rc,spec_ssa_eof,spec_nlsa_eof,gaussian_ssa,gaussian_nlsa,li_harmonics_ssa,li_harmonics_nlsa,ssa_trend_harm,nlsa_trend_harm,freq_ssa,freq_nlsa,ssa_harm_var,nlsa_harm_var,spec_ssa,spec_res_ssa,spec_nlsa,spec_res_nlsa = p

    modenumber = 12
    modenumber_k = 48
    xstep_k = 16
    k_xticks = (xstep_k:xstep_k:modenumber_k,string.(xstep_k:xstep_k:modenumber_k))
    spec_xticks = (1:6,string.(1:6))
    smallfs = fontsize-6
    fs = fontsize
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
    xticklabelsize = fs-2,
    yticklabelsize = fs-2,
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
    xticklabelsize = fs-2,
    yticklabelsize = fs-2,
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
    xticklabelsize = fs-2,
    yticklabelsize = fs-2,
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
    xticklabelsize = fs-2,
    yticklabelsize = fs-2,
    ylabelsize = fs,)

    hideydecorations!(ax_modes,ticks = false,ticklabels=false,grid=false)

    ax_freq = Axis(F[5:9,3:4],
    limits=freqlimits,#,yscale=log10
    xticks=1:7,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency [1/a]",
    ylabel="relative power",
    xlabelsize = fs,
    xticklabelsize = fs-2,
    yticklabelsize = fs-2,
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

function large_mode_figure_flags(savedirname,data_tensor,outdir)
    F = Figure(resolution=(2400,1000))

    ga = F[1:12, 1:8] = GridLayout()
    gb = F[1:12, 9:16] = GridLayout()
    gc = F[1:12, 17:24] = GridLayout()
    
    #gl1 = F[1,0] = GridLayout()
    #gl2 = F[2,0] = GridLayout()
    #gl3 = F[7,0] = GridLayout()

    gt1 = F[0,2] = GridLayout()
    gt2 = F[0,10] = GridLayout()
    gt3 = F[0,18] = GridLayout()

    #gl = F[8:12,0] = GridLayout()

    """
    local, filtered coordinates!
    """
    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]
    #measurement
    spot = 9
    var = 3
    p = local_parameters(spots[spot],vars[var],outdir)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(ga,p,rich("GPP [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor)
    title1 = "(a) $(variables_names[vars[var]]) $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])"

    #resolved
    spot = 2
    var = 3
    p = local_parameters(spots[spot],vars[var],outdir)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(gb,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor)
    title2 = "(b) $(variables_names[vars[var]]) $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])"

    #unresolved
    spot = 1
    var = 2
    p = local_parameters(spots[spot],vars[var],outdir)
    flag = flags[:,spots[spot],[6,7,5,2,9,8][var]]
    mode_figure_flags(gc,p,rich("RECO [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor)
    title3 = "(c) $(variables_names[vars[var]]) $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])"

    for (label, layout) in zip([title1,title2,title3], [gt1,gt2,gt3])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font = :bold,
            #padding = (0, 5, 5, 0),
            #halign = :right
            )
    end

    """
    bigfs = 26
    label1 = "signal and\n reconstruction"
    label2 = "relative power"
    label3 = "individual modes"
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
    """

    for g in [ga, gb, gc]
        g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end

    for g in [gt1,gt2,gt3]
        colgap!(g, 0)
        rowgap!(g, 0)
        g.alignmode = Mixed(right = 0,top=0,bottom=0,left=0)
    end

    #colsize!(F.layout, 0, 15)
    #colsize!(F.layout, 2, 75)
    #colsize!(F.layout, 3, 75)
    #colsize!(F.layout, 4, 75)

    save(savedirname,F)
end

#heatmap
function combined_heatmap(savedirname)

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
    fs = fontsize
    fs_plot = fontsize-4


    function plot_heatmap(g,data)

        data = data[spots,vars]

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
                xticks = (1:length(spots),spotslist[spots].*" ".*IGBP_reduced),
                xticklabelrotation = pi/2,
                #xlabel="site",
                #ylabel = "variable",
                xlabelsize = fs,
                xticklabelsize = fs-4,
                yticks = (1:length(varnames),varnames),
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

    ax_a,hm_a = plot_heatmap(g_a,ssa_h_raw)
    ax_b,hm_b = plot_heatmap(g_b,nlsa_h_raw)
    ax_c,hm_c = plot_heatmap(g_c,ssa_h_6)
    ax_d,hm_d = plot_heatmap(g_d,nlsa_h_6)
    ax_e,hm_e = plot_heatmap(g_e,ssa_h_4)
    ax_f,hm_f = plot_heatmap(g_f,nlsa_h_4)
    ax_g,hm_g = plot_heatmap(g_g,ssa_h_3)
    ax_h,hm_h = plot_heatmap(g_h,nlsa_h_3)

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

    save(savedirname,F)
end