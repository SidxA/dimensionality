# Function to create a motivation plot
function motivation_plot(savedirname, outdir)
    # Set spot and variable indices
    spot = 9
    vari = 9

    # Local parameters calculation
    p = local_parameters(spot, vari, data_tensor, outdir)

    # Function to compute individual harmonics
    function individual_harmonics(p, rescale = false)
        # Retrieve parameters
        spot, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p
        
        # Calculate SSA and NLSA trend reconstructed components
        l = embed_lag(signal, W)'
        ssa_trend_rc = hcat([reconstructor((l * ssa_Eof)[:, i], ssa_Eof[:, i], N, W) for i in li_harmonics_ssa]...)
        nlsa_trend_rc = hcat([reconstructor((l * nlsa_Eof)[:, i], nlsa_Eof[:, i], N, W) for i in li_harmonics_nlsa]...)
        
        # Rescale if specified
        if rescale == true
            m = means[spot, vari]
            s = stds[spot, vari]
            ssa_trend_rc = hcat([back_trafo(ssa_trend_rc[:, i], m, s) for i in 1:size(ssa_trend_rc, 2)]...)
            nlsa_trend_rc = hcat([back_trafo(nlsa_trend_rc[:, i], m, s) for i in 1:size(nlsa_trend_rc, 2)]...)
        end
        
        return ssa_trend_rc, nlsa_trend_rc
    end

    # Back transformation function
    back_trafo(data, mean, std) = (data .* std) .+ mean

    # Compute individual harmonics
    ssa_trend_rc, nlsa_trend_rc = individual_harmonics(p, false)

    # Rescale local parameters
    p = rescale_local_parameters(p, data_tensor)

    # Retrieve parameters
    spot, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p
    
    # Compute mean and standard deviation of raw data
    m = mean(data_raw[:, spot, vari])
    s = std(data_raw[:, spot, vari])

    # Compute fundamental components
    ssa_fund = back_trafo(sum(ssa_trend_rc[:, 1:2], dims=2)[:], m, s)
    nlsa_fund = back_trafo(sum(nlsa_trend_rc[:, 1:2], dims=2)[:], m, s)

    # Create string for seasonal cycle information
    seasonalitystring = ""
    seasonalitystring *= "ssa f \t" * string(round.(freq_ssa, digits=1)) * "\n"
    seasonalitystring *= "nlsa f \t" * string(round.(freq_nlsa, digits=1)) * "\n"
    seasonalitystring *= "ssa var \t" * string(round.(ssa_harm_var, digits=1)) * "\n"
    seasonalitystring *= "nlsa var \t" * string(round.(nlsa_harm_var, digits=1)) * "\n"

    # Create figure
    offset = 15
    tiks = Array(-2:2:2)
    yticks = (cat(tiks, tiks .+ offset, dims=1), cat(string.(tiks), string.(tiks), dims=1))
    F = Figure(resolution=(800, 400))
    ax_time = Axis(F[1, 1], xticks=Int.(floor.(years[1]:3:years[end])), yticks=yticks, yminorgridvisible=true, yminorgridstyle=:dot, yminorticks=IntervalsBetween(5), xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(3), xlabel="time (a)", ylabel=rich("NEE [gC m", superscript("-2"), "d", superscript("-1"), "]"), xlabelsize=fontsize, ylabelsize=fontsize, xticklabelsize=fontsize-4, yticklabelsize=fontsize-4)

    # Signal background
    scatter!(ax_time, years, signal .+ offset, linewidth=lw, color=color_signal, marker=:x, markersize=ms)
    li = lines!(ax_time, years, signal .+ offset, linewidth=lw_s, color=color_signal, label="signal")

    # SSA trend
    li = lines!(ax_time, years, ssa_trend_harm .+ offset, linewidth=lw, color=color_ssa, linestyle=:solid, label="SSA")

    # NLSA trend
    li = lines!(ax_time, years, nlsa_trend_harm .+ offset, linewidth=lw, color=color_nlsa, linestyle=:solid, label="NLSA")

    # Signal fundamental
    scatter!(ax_time, years, signal, linewidth=lw, color=color_signal, marker=:x, markersize=ms)
    lines!(ax_time, years, signal, linewidth=lw_s, color=color_signal)

    # SSA fundamental
    lines!(ax_time, years, ssa_fund, linewidth=lw, color=color_ssa, linestyle=:solid)

    # NLSA fundamental
    lines!(ax_time, years, nlsa_fund, linewidth=lw, color=color_nlsa, linestyle=:solid)

    # Add text annotations
    text!(ax_time, 0, 1, text="seasonal cycle", align=(:left, :top), offset=(4, -2), space=:relative, fontsize=fontsize)
    text!(ax_time, 0, 0.02, text="only fundamental", align=(:left, :bottom), offset=(4, -2), space=:relative, fontsize=fontsize)
    text!(ax_time, 0.95, 1, text="$(spotslist[spot]) $(igbpclass) $(varname)", align=(:right, :top), space=:relative, fontsize=fontsize)

    # Customize axis and legend
    hidespines!(ax_time, :t, :r)
    axislegend(ax_time, tellheight=false, tellwidth=false, margin=(4, 4, 4, 4), halign=:right, valign=:bottom, orientation=:horizontal, labelsize=fontsize-2)

    # Save the figure
    save(savedirname, F)
end
