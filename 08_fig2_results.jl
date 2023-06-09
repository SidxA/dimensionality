include("/net/home/lschulz/dimensionality/06_load_figuredata.jl")

"""
the following figures are built using the reduced number of spots and variables
"""

# choose an unfiltered time series
data_tensor = data_raw
outdir = outdir_raw

#choose a filtered time series
#data_tensor = data_filtered
#outdir = outdir_filtered

#data saving directory
dir = init_logging()
savedirname = dir * "test.png"


#individual overview plot showing the time series, the spectrum, the modes and the spectra of the modes
#need to be given a Figure() as F in order to combine them to a single figure
#varname resolved is the name of the variable together with the corresponding unit, making up the time series axis label
function mode_figure_flags(F,p,varname_resolved,flags,data_tensor,ylim)

    # Rescale local parameters using data tensor
    p = rescale_local_parameters(p, data_tensor)

    # Extract variables from parameter tuple
    spot, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p

    # Define constants
    modenumber = 12  # Total number of modes
    modenumber_k = 16  # Total number of k modes
    xstep_k = 4  # Step size for k modes
    k_xticks = (xstep_k:xstep_k:modenumber_k, string.(xstep_k:xstep_k:modenumber_k))  # X-axis ticks for k modes
    spec_xticks = (1:7, string.(1:7))  # X-axis ticks for spectral analysis
    smallfs = fontsize - 6  # Small font size
    fs = fontsize  # Font size
    speclimits = (freq_domain_N[1], 7.05, 10^-4, 4)  # Limits for spectral analysis plot
    klimits = (0, modenumber_k + 1, 10^-4, 4)  # Limits for k modes plot
    modeslimits = (0, 2556 / 365.25, -0.04, (2 * modenumber + 2.5) * 0.08)  # Limits for modes plot
    modeyticks = (vcat(0.08:0.08*2:(modenumber)*0.08, (modenumber + 2)*0.08:0.08*2:(2*modenumber + 1)*0.08), string.(vcat(1:2:modenumber, 1:2:modenumber)))  # Y-axis ticks for modes plot
    freqlimits = (1/10, 7, -0.2, modenumber*2 + 2.5)  # Limits for frequency plot
    gapsize = 8  # Size of gap
    lw = 3  # Line width
    ms = 12  # Marker size
    scatterstep = 10  # Step size for scatter plot
    color_signal = "darkorange"  # Color for signal

    #figure
    # Create an axis for time plot
    ax_time = Axis(F[1:2,1:4],
        xticks=Int.(floor.(years[1]:3:years[end])),  # Set x-axis ticks as integers based on the floor of years in intervals of 3
        limits=(years[1], years[end], ylim[1],ylim[2]),  # Set the limits of the axis based on years and signal values
        xminorticksvisible=true,  # Display minor ticks on the x-axis
        xminorgridvisible=true,  # Display minor grid lines on the x-axis
        xminorticks=IntervalsBetween(3),  # Set the interval between minor ticks as 3
        xlabel="time [a]",  # Set x-axis label
        ylabel=varname_resolved,  # Set y-axis label as varname_resolved
        xlabelsize=fs,  # Set x-axis label font size
        xticklabelsize=fs-2,  # Set x-axis tick label font size
        yticklabelsize=fs-2,  # Set y-axis tick label font size
        ylabelsize=fs,)  # Set y-axis label font size
        #alignmode = Outside())

    hideydecorations!(ax_time, ticks=false, ticklabels=false, grid=false, label=false)  # Hide y-axis decorations (ticks, tick labels, grid, and label)

    # Create an axis for spectral analysis plot
    ax_spec = Axis(F[3:4,1:3], yscale=log10,
        limits=speclimits,  # Set the limits of the axis based on speclimits
        xticks=spec_xticks,  # Set x-axis ticks based on spec_xticks
        xminorticksvisible=true,  # Display minor ticks on the x-axis
        xminorgridvisible=true,  # Display minor grid lines on the x-axis
        xminorticks=IntervalsBetween(4),  # Set the interval between minor ticks as 4
        xlabel="frequency [1/a]",  # Set x-axis label
        ylabel="relative power",  # Set y-axis label
        xlabelsize=fs,  # Set x-axis label font size
        xticklabelsize=fs-2,  # Set x-axis tick label font size
        yticklabelsize=fs-2,  # Set y-axis tick label font size
        ylabelsize=fs,)  # Set y-axis label font size
        #alignmode = Outside())

    hideydecorations!(ax_spec, ticks=false, ticklabels=false, grid=false)  # Hide y-axis decorations (ticks, tick labels, grid)

    # Create an axis for k modes plot
    ax_k = Axis(F[3:4,4], yscale=log10,
        limits=klimits,  # Set the limits of the axis based on klimits
        xminorticksvisible=true,  # Display minor ticks on the x-axis
        xminorgridvisible=true,  # Display minor grid lines on the x-axis
        xminorticks=IntervalsBetween(2),  # Set the interval between minor ticks as 2
        xlabel="dim",  # Set x-axis label
        xlabelsize=fs,  # Set x-axis label font size
        xticks=k_xticks,  # Set x-axis ticks based on k_xticks
        xticklabelsize=fs-2,  # Set x-axis tick label font size
        yticklabelsize=fs-2,  # Set y-axis tick label font size
        ylabelsize=fs,)  # Set y-axis label font size
        #alignmode = Outside())

    hideydecorations!(ax_k, grid=false)  # Hide y-axis decorations (grid)

    # Create an axis for individual modes plot
    ax_modes = Axis(F[5:9,1:2],
        xticks=0:Int(floor(W/365.25)),  # Set x-axis ticks from 0 to the floor of W divided by 365.25
        yticks=modeyticks,  # Set y-axis ticks based on modeyticks
        limits=modeslimits,  # Set the limits of the axis based on modeslimits
        xminorticksvisible=true,  # Display minor ticks on the x-axis
        xminorgridvisible=true,  # Display minor grid lines on the x-axis
        xminorticks=IntervalsBetween(2),  # Set the interval between minor ticks as 2
        xlabel="time [a]",  # Set x-axis label
        ylabel="individual modes",  # Set y-axis label
        xlabelsize=fs,  # Set x-axis label font size
        xticklabelsize=fs-2,  # Set x-axis tick label font size
        yticklabelsize=fs-2,  # Set y-axis tick label font size
        ylabelsize=fs,)  # Set y-axis label font size
        #alignmode = Outside())

    hideydecorations!(ax_modes, ticks=false, ticklabels=false, grid=false)  # Hide y-axis decorations (ticks, tick labels, grid)

    # Create an axis for frequency plot
    ax_freq = Axis(F[5:9,3:4],
        limits=freqlimits,  # Set the limits of the axis based on freqlimits
        xticks=1:7,  # Set x-axis ticks from 1 to 7
        yticksvisible=false,  # Hide y-axis ticks
        yticklabelsvisible=false,  # Hide y-axis tick labels
        xlabel="frequency [1/a]",  # Set x-axis label
        ylabel="relative power",  # Set y-axis label
        xlabelsize=fs,  # Set x-axis label font size
        xticklabelsize=fs-2,  # Set x-axis tick label font size
        yticklabelsize=fs-2,  # Set y-axis tick label font size
        ylabelsize=fs,)  # Set y-axis label font size
        #alignmode = Outside())

    hideydecorations!(ax_freq)  # Hide y-axis decorations




    #flags
    vlines!(ax_time,years,color=flags,ymin = 0.0, ymax = 1.0,colormap = :grays)
    
    #plotting
    #scatter!(ax_time,years[1:scatterstep:end],signal[1:scatterstep:end],color=color_signal,markersize = ms,marker=:x)
    # Add lines to the ax_time axis
    signal_l = lines!(ax_time, years, signal, color=color_signal, linewidth=lw, label="signal")
    ssa_l = lines!(ax_time, years, ssa_trend_harm, color=color_ssa, linewidth=lw, label="SSA")
    nlsa_l = lines!(ax_time, years, nlsa_trend_harm, color=color_nlsa, linewidth=lw, label="NLSA")

    # Add lines to the ax_spec axis
    lines!(ax_spec, freq_domain_N, spec_signal, color=color_signal, linewidth=lw, label="signal")
    signal_s = scatter!(ax_spec, freq_domain_N, spec_signal, color=color_signal, markersize=ms, marker=:x)
    #lines!(ax_spec, freq_domain_N, spec_ssa_rc, color=color_ssa, linewidth=lw, label="SSA")
    #lines!(ax_spec, freq_domain_N, spec_nlsa_rc, color=color_nlsa, linewidth=lw, label="NLSA")
    lines!(ax_spec, freq_domain_N, spec_ssa, color=color_ssa, linewidth=lw, label="SSA")
    lines!(ax_spec, freq_domain_N, spec_nlsa, color=color_nlsa, linewidth=lw, label="NLSA")
    
    # Add lines to the ax_k axis
    lines!(ax_k, 1:modenumber_k, ssa_cap_var[1:modenumber_k], color=color_ssa, linewidth=lw)
    lines!(ax_k, 1:modenumber_k, nlsa_cap_var[1:modenumber_k], color=color_nlsa, linewidth=lw)
    #add markers at harmonics
    scatter!(ax_k, li_harmonics_ssa, ssa_cap_var[li_harmonics_ssa], color=color_ssa, markersize=ms+2, marker=:x)
    scatter!(ax_k, li_harmonics_nlsa, nlsa_cap_var[li_harmonics_nlsa], color=color_nlsa, markersize=ms+2, marker=:x)




    #modes 
    # Define color for fit lines
    color_fit = "grey"

    # Iterate over the modes
    for i = 1:modenumber
        # Get mode and Fourier transform values for SSA
        mode = ssa_Eof[:,i]
        Four = spec_ssa_eof[:,i]
        years_modes = (1:length(mode)) ./ 365.25

        # Add lines to ax_modes and ax_freq axes for SSA mode
        lines!(ax_modes, years_modes, mode .+ (i * 0.08), color=color_ssa)
        lines!(ax_freq, freqs_w, abs.(Four) .+ i, color=color_ssa)

        # Get mode and Fourier transform values for NLSA
        mode = nlsa_Eof[:,i]
        Four = spec_nlsa_eof[:,i]

        # Add lines to ax_modes and ax_freq axes for NLSA mode
        lines!(ax_modes, years_modes, mode .+ ((modenumber + 1 + i) * 0.08), color=color_nlsa)
        lines!(ax_freq, freqs_w, abs.(Four) .+ i .+ modenumber .+ 1, color=color_nlsa)
    end

    # Generate boxstring for displaying height and frequency information
    boxstring = "f \t height \n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_nlsa[k][1], digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_nlsa[k][2], digits=1))
        boxstring *= "\n"
    end
    boxstring *= "\n"
    for k = modenumber:-1:1
        boxstring *= string(round(gaussian_ssa[k][1], digits=1))
        boxstring *= "\t"
        boxstring *= string(round(gaussian_ssa[k][2], digits=1))
        boxstring *= "\n"
    end

    # Link y-axes of ax_spec and ax_k axes
    linkyaxes!(ax_spec, ax_k)


    return [ax_time,ax_spec,ax_k,ax_modes,ax_freq]
end




function examples_regular(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 8
    var = 1
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [-1,20]
    axes = mode_figure_flags(ga,p,rich("GPP [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("GPP [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "GPP $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

function examples_resolved(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 2
    var = 3
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [-8,3]
    axes = mode_figure_flags(ga,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "NEE $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

function examples_unresolved(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 5
    var = 3
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [-5,5]
    axes = mode_figure_flags(ga,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("NEE [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "NEE $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

function examples_deficient(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 7
    var = 2
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [-2,15]
    axes = mode_figure_flags(ga,p,rich("RECO [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("RECO [gC m", superscript("-2")," d",superscript("-1"),"]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "RECO $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

function examples_appendix1(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 5
    var = 4
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [0,400]
    axes = mode_figure_flags(ga,p,rich("SW_IN [W m", superscript("-2"),"]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("SW_IN [W m", superscript("-2"),"]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "SW_IN $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

function examples_appendix2(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 7
    var = 6
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [0,30]
    axes = mode_figure_flags(ga,p,rich("TS [deg C]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("TS [deg C]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "TS $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

function examples_appendix3(savedirname)

    F = Figure(resolution=(1600,1000))
    colsize = 655
    rowsize1 = 50
    rowsize2 = 840

    gt1 = F[1,2] = GridLayout(1,1, alignmode = Outside())
    gt2 = F[1,3] = GridLayout(1,1, alignmode = Outside())

    ga = F[2, 2] = GridLayout(1,1,alignmode = Inside())
    gb = F[2, 3] = GridLayout(1,1,alignmode = Inside())
    
    gl1 = F[1,1] = GridLayout(1,1,alignmode = Outside())
    gl2 = F[2,1] = GridLayout(1,1,alignmode = Outside())


    spots = mask_IGBP(IGBP_list)[1]
    vars = mask_vari(variables_names)[1]

    #unfiltered
    spot = 2
    var = 7
    p = local_parameters(spots[spot],vars[var],outdir_raw)
    flag = flags[:,spots[spot],[6,7,5,2,1,9,8][var]]
    ylim = [0,50]
    axes = mode_figure_flags(ga,p,rich("SWC [%]"),flag,data_tensor,ylim)
    title1 = "(a) unfiltered  "

    #filtered
    p = local_parameters(spots[spot],vars[var],outdir_f6)
    axes = mode_figure_flags(gb,p,rich("SWC [%]"),flag,data_tensor,ylim)
    title2 = "(b) filtered  "

    [hideydecorations!(ax,grid=:false,ticks=:false) for ax in axes]

    colsize!(F.layout, 1, 100)
    colsize!(F.layout, 2, colsize)
    colsize!(F.layout, 3, colsize)
    rowsize!(F.layout,1,rowsize1)
    rowsize!(F.layout,2,rowsize2)


    rowgap!(F.layout,0)
    colgap!(F.layout,1,5)
    colgap!(F.layout,2,20)

    Label(gl2[1,2], "modes            spectra          signal  ",
    rotation=pi/2,fontsize=fontsize+4,font=:bold,halign=:right)
    #Label(gl2[2,2], "spectra",rotation=pi/2,fontsize=fontsize+4,font=:bold)
    #Label(gl2[3,2], "      modes",rotation=pi/2,fontsize=fontsize+4,font=:bold)

    Label(gl2[1,1],
    "SW_IN $(spotslist[spots[spot]]) $(IGBP_list[spots[spot]])",
    rotation=pi/2,
    font=:bold,
    fontsize=fontsize+4,
    halign=:right)

    elem_1 = [LineElement(color = "darkorange", linestyle = :solid,linewidth = 10)]
    elem_2 = [LineElement(color = :darkgreen, linestyle = :solid,linewidth = 10)]
    elem_3 = [LineElement(color = :purple, linestyle = :solid,linewidth = 10)]
    elem_4 = MarkerElement(color = [:black,:grey90,:grey80], marker = :vline, markersize = 40,
        points = Point2f[(0.0, 0.4), (0.4, 0.4), (0.8, 0.4)])

    Legend(gl2[7,1:2],
    [elem_1, elem_2, elem_3, elem_4],
    ["Signal", "SSA", "NLSA", "Q flag"],
    labelsize = fontsize+4,padding=(20,20,0,0))

    Label(gl2[8,1], "   ",rotation=pi/2,fontsize=fontsize+4)


    for (label, layout) in zip([title1,title2], [gt1,gt2])
        Label(layout[1, 1], label,
            fontsize = fontsize +4,
            font=:bold,
            padding = (0, 5, 20, 20),
            halign = :right,
            justification = :right
            )
    end

    save(savedirname,F)

end

#examples_regular(dir*"results/regular.pdf")
#examples_resolved(dir*"results/resolved.pdf")
#examples_unresolved(dir*"results/unresolved.pdf")
#examples_deficient(dir*"results/deficient.pdf")
#examples_appendix1(dir*"results/appendix1.pdf")
#examples_appendix2(dir*"results/appendix2.pdf")
#examples_appendix3(dir*"results/appendix3.pdf")
