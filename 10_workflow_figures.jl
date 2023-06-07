include("/net/home/lschulz/dimensionality/06_load_figuredata.jl")

# choose an unfiltered time series
data_tensor = data_raw
outdir = outdir_raw

#data saving directory
dir = init_logging()
savedirname = dir * "test.png"

color_signal = "black"
color_background = "gray"
# Set spot and variable indices
spot = 9
vari = 9

# Local parameters calculation
p = local_parameters(spot, vari, outdir)

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

# Rescale local parameters
p = rescale_local_parameters(p, data_tensor)

# Compute individual harmonics
ssa_trend_rc, nlsa_trend_rc = individual_harmonics(p, false)

# Retrieve parameters
spot, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p

"""
# Create figure for time series
"""
F = Figure(resolution=(800, 400))
ax_time = Axis(F[1, 1], 
xticks=Int.(floor.(years[1]:3:years[end])), 
yminorgridvisible=true, 
yminorgridstyle=:dot, 
yminorticks=IntervalsBetween(5), 
xminorticksvisible=true, 
xminorgridvisible=true, 
xminorticks=IntervalsBetween(3), 
xlabel="time [a]", 
ylabel=rich("NEE [gC m", superscript("-2"), "d", superscript("-1"), "]"), 
xlabelsize=fontsize, 
ylabelsize=fontsize, 
xticklabelsize=fontsize-4, 
yticklabelsize=fontsize-4)

# Signal background
scatter!(ax_time, years, signal, linewidth=lw, color=color_signal, marker=:x, markersize=ms)
lines!(ax_time, years, signal, linewidth=lw_s, color=color_signal, label="signal")

# Customize axis and legend
hidespines!(ax_time, :t, :r)
# Save the figure
save(dir*"workflow/time_series.pdf", F)
hideydecorations!(ax_time)
hidespines!(ax_time, :l)
ax_time.xticks = 2008:4:2020
save(dir*"workflow/time_series_minimal.pdf", F)
hidespines!(ax_time)
hidedecorations!(ax_time)
save(dir*"workflow/time_series_nospines.pdf", F)


"""
# Create figure for seasonal cycle
"""
F = Figure(resolution=(800, 400))
ax_time = Axis(F[1, 1], 
xticks=Int.(floor.(years[1]:3:years[end])), 
yminorgridvisible=true, 
yminorgridstyle=:dot, 
yminorticks=IntervalsBetween(5), 
xminorticksvisible=true, 
xminorgridvisible=true, 
xminorticks=IntervalsBetween(3), 
xlabel="time [a]", 
ylabel=rich("NEE [gC m", superscript("-2"), "d", superscript("-1"), "]"), 
xlabelsize=fontsize, 
ylabelsize=fontsize, 
xticklabelsize=fontsize-4, 
yticklabelsize=fontsize-4)

# Signal background
scatter!(ax_time, years, signal, linewidth=lw, color=color_background, marker=:x, markersize=ms)
lines!(ax_time, years, signal, linewidth=lw_s, color=color_background, label="signal")

lines!(ax_time,years,sum(nlsa_trend_rc, dims=2)[:], linewidth=lw, color=color_signal, label="SSA trend")
# Customize axis and legend
hidespines!(ax_time, :t, :r)
# Save the figure
save(dir*"workflow/seasonal_cycle.pdf", F)
hideydecorations!(ax_time)
hidespines!(ax_time, :l)
ax_time.xticks = 2008:4:2020
save(dir*"workflow/seasonal_cycle_minimal.pdf", F)
hidespines!(ax_time)
hidedecorations!(ax_time)
save(dir*"workflow/seasonal_cycle_nospines.pdf", F)

"""
# Create figure for W = 7 a
"""
time_interval = 1:7*365
F = Figure(resolution=(800, 400))
ax_time = Axis(F[1, 1], 
#xticks=Int.(floor.(years[1]:3:years[end])), 
yminorgridvisible=true, 
yminorgridstyle=:dot, 
yminorticks=IntervalsBetween(5), 
xminorticksvisible=true, 
xminorgridvisible=true, 
xminorticks=IntervalsBetween(3), 
xlabel="time [a]", 
ylabel=rich("NEE [gC m", superscript("-2"), "d", superscript("-1"), "]"), 
xlabelsize=fontsize, 
ylabelsize=fontsize, 
xticklabelsize=fontsize-4, 
yticklabelsize=fontsize-4)

lines!(ax_time,(time_interval)./365.25,signal[time_interval], linewidth=lw, color=color_signal, label="SSA trend")
# Customize axis and legend
hidespines!(ax_time, :t, :r)
# Save the figure
save(dir*"workflow/embedding.pdf", F)
hideydecorations!(ax_time)
hidespines!(ax_time, :l)
ax_time.xticks = 1:3:7
save(dir*"workflow/embedding_minimal.pdf", F)
hidespines!(ax_time)
hidedecorations!(ax_time)
save(dir*"workflow/embedding_nospines.pdf", F)

"""
# Create figure for 7 a of sine
"""
time_interval = 1:7*365
F = Figure(resolution=(800, 400))
ax_time = Axis(F[1, 1], 
#xticks=Int.(floor.(years[1]:3:years[end])), 
yminorgridvisible=true, 
yminorgridstyle=:dot, 
yminorticks=IntervalsBetween(5), 
xminorticksvisible=true, 
xminorgridvisible=true, 
xminorticks=IntervalsBetween(3), 
xlabel="time [a]", 
ylabel=rich("NEE [gC m", superscript("-2"), "d", superscript("-1"), "]"), 
xlabelsize=fontsize, 
ylabelsize=fontsize, 
xticklabelsize=fontsize-4, 
yticklabelsize=fontsize-4)

lines!(ax_time,(time_interval)./365.25,sin.(2* pi*time_interval ./ 365), linewidth=lw, color=color_signal, label="SSA trend")
# Customize axis and legend
hidespines!(ax_time, :t, :r)
# Save the figure
save(dir*"workflow/fundamental.pdf", F)
hideydecorations!(ax_time)
hidespines!(ax_time, :l)
ax_time.xticks = 1:3:7
save(dir*"workflow/fundamental_minimal.pdf", F)
hidespines!(ax_time)
hidedecorations!(ax_time)
save(dir*"workflow/fundamental_nospines.pdf", F)

"""
modes no
"""
data_tensor = data_f4
outdir = outdir_f4

spot = spots[2]
vari = vars[3]

# Local parameters calculation
p = local_parameters(spot, vari, outdir)

# Retrieve parameters
spot, W, vari, years, varname, igbpclass, freq_domain_N, freq_domain_w, freqs_w, freqs, signal, ssa_Eof, nlsa_Eof, nlsa_eps, ssa_rec, nlsa_rec, ssa_cap_var, nlsa_cap_var, spec_signal, spec_ssa_rc, spec_nlsa_rc, spec_ssa_eof, spec_nlsa_eof, gaussian_ssa, gaussian_nlsa, li_harmonics_ssa, li_harmonics_nlsa, ssa_trend_harm, nlsa_trend_harm, freq_ssa, freq_nlsa, ssa_harm_var, nlsa_harm_var, spec_ssa, spec_res_ssa, spec_nlsa, spec_res_nlsa = p

F = Figure()
ax = Axis(F[1,1],limits = (0.5,3.5,0.01,1),
xlabel=L"f \text{ [a}^{-1}]", 
ylabel="relative power", 
xlabelsize=fontsize, 
ylabelsize=fontsize, 
xticklabelsize=fontsize-4, 
yticklabelsize=fontsize-4,
xticks = [1,2,3],
xtickwidth = 2,
xticksize = 10,)
lines!(ax,freqs_w,spec_ssa_eof[:,5],color = "black",linewidth = 6)

# Customize axis and legend
hidespines!(ax, :t, :r)
# Save the figure
save(dir*"workflow/modes_no.pdf", F)
hideydecorations!(ax)
hidespines!(ax, :l)
ax.xticks = 1:1:7
save(dir*"workflow/modes_no_minimal.pdf", F)
hidespines!(ax)
hidedecorations!(ax)
save(dir*"workflow/modes_no_nospines.pdf", F)

"""
modes yes
"""

F = Figure()
ax = Axis(F[1,1],limits = (0.5,3.5,0.01,1),
xlabel=L"f \text{ [a}^{-1}]", 
ylabel="relative power", 
xlabelsize=fontsize, 
ylabelsize=fontsize, 
xticklabelsize=fontsize-4, 
yticklabelsize=fontsize-4,
xticks = [1,2,3],
xtickwidth = 2,
xticksize = 10,)
lines!(ax,freqs_w,spec_nlsa_eof[:,5],color = "black",linewidth = 6)

# Customize axis and legend
hidespines!(ax, :t, :r)
# Save the figure
save(dir*"workflow/modes_yes.pdf", F)
hideydecorations!(ax)
hidespines!(ax, :l)
ax.xticks = 1:1:7
save(dir*"workflow/modes_yes_minimal.pdf", F)
hidespines!(ax)
hidedecorations!(ax)
save(dir*"workflow/modes_yes_nospines.pdf", F)

