"""
the following figures are built using the reduced number of spots and variables
"""


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
    spot = 8
    var = 1
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

