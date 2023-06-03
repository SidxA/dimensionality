include("/net/home/lschulz/dimensionality/06_load_figuredata.jl")


#data saving directory
dir = init_logging()
savedirname = dir * "test.png"


#characterization map
function characterization_heatmap(savedirname)

    function plot_heatmap(g, data)
        # Subset the data based on spots and variables
        data = data[spots, vars]
    
        # Replace values less than 2 with NaN
        data[data .< 2] .= NaN
    
        x = 1:size(data, 1)
        y = 1:size(data, 2)
    
        # Create a heatmap plot on the specified grid position
        ax, hm = heatmap(g[1, 1], x, y, data,
            levels = 2:1:highclip,
            colormap = cgrad(:heat,highclip-1,categorical = true),
            colorrange = (1.5, highclip+0.5),
            framevisible = false,
            axis = (
                xticks = (1:length(spots), spotslist[spots] .* " " .* IGBP_reduced),
                xticklabelrotation = pi/2,
                xlabelsize = fs,
                xticklabelsize = fs - 4,
                yticks = (1:length(varnames), varnames),
                yticklabelsize = fs - 4,
                ylabelsize = fs,
            ))
    
        hidespines!(ax, :t, :r)  # Hide top and right spines
    
        return ax, hm
    end
    


    function write_symbol(ax, data, symbol, x_shift, y_shift) # 5bins
        symbol_fontsize_list = [0, 8, 12, 16, 20]  # Font size options for different symbol levels
        for i in 1:size(data)[1]
            str = data[i, :]
            map(enumerate(str)) do (j, c)
                fontsize = symbol_fontsize_list[Int(c)]  # Determine the font size based on the symbol level
                text!(ax, i + x_shift, j + y_shift, text=symbol,
                    align = (:center, :center),
                    color = :black,
                    fontsize = fontsize,
                )
            end
        end
    end
    
    function artifact_symbol(ax, data, x_shift, y_shift) #boolean
        symbol = "X"  # Symbol to represent artifact
        fs = 20  # Font size for the symbol
        for i in 1:size(data)[1]
            str = data[i, :]
            map(enumerate(str)) do (j, c)
                if c  # If the artifact is present
                    text!(ax, i + x_shift, j + y_shift, text=symbol,
                        align = (:center, :center),
                        color = :black,
                        fontsize = fs,
                    )
                end
            end
        end
    end
    

    function scatter_symbol(ax, data, symbol, x_shift, y_shift)
        symbol_sizes = [0, 8, 16]  # Sizes for the scatter symbols
        for i in 1:size(data)[1]
            str = data[i, :]
            map(enumerate(str)) do (j, c)
                scatter!(ax, i + x_shift, j + y_shift, markersize=symbol_sizes[Int(c)], marker=symbol,
                    align = (:center, :center),
                    color = :black,
                )
            end
        end
    end
    



    F = Figure(resolution = (600,900))
    
    g_a = F[1:3,1:10] = GridLayout()
    g_b = F[1:3,11:20] = GridLayout()
    g_c = F[4:6,1:10] = GridLayout()
    g_d = F[4:6,11:20] = GridLayout()
    g_e = F[7:9,1:10] = GridLayout()
    g_f = F[7:9,11:20] = GridLayout()

    g_cbar = F[0,1:20] = GridLayout()

    g_lab_1 = F[2,0] = GridLayout()
    g_lab_2 = F[5,0] = GridLayout()
    g_lab_3 = F[8,0] = GridLayout()


    highclip = 8
    fs = fontsize
    fs_plot = fontsize-4

    for g in [g_a,g_b,g_c,g_d,g_e,g_f]
        colgap!(g, 0)
        rowgap!(g, 0)
        g.alignmode = Mixed(right = 0,top=0,bottom=0)
    end

    ax_a,hm_a = plot_heatmap(g_a,ssa_h_raw)
    ax_b,hm_b = plot_heatmap(g_b,nlsa_h_raw)
    ax_c,hm_c = plot_heatmap(g_c,ssa_h_6)
    ax_d,hm_d = plot_heatmap(g_d,nlsa_h_6)
    ax_e,hm_e = plot_heatmap(g_e,ssa_h_4)
    ax_f,hm_f = plot_heatmap(g_f,nlsa_h_4)


    #on the right, top2
    for (ax,g) in zip([ax_b,ax_d],[g_b,g_d])
        hidedecorations!(ax,grid = false,ticks=false)
    end

    #on the left, top 2
    for (ax,g) in zip([ax_a,ax_c],[g_a,g_c])
        hidexdecorations!(ax,grid = false,ticks=false)
    end

    #bottom figures
    hideydecorations!(ax_f,grid = false,ticks=false)


    Colorbar(g_cbar[1,1:6], hm_a, vertical = false,ticklabelsize=fontsize,
    )

    subfigure_labels = [
        L"\text{(a) SSA      }",
        L"\text{(b) NLSA      }",
        L"\text{(c) SSA           } f \leq 6/\text{a}",
        L"\text{(d) NLSA        } f \leq 6/\text{a}",
        L"\text{(e) SSA           } f \leq 4/\text{a}",
        L"\text{(f) NLSA         } f \leq 4/\text{a}"
    ]
    gs = [g_a,g_b,g_c,g_d,g_e,g_f]

    
    for (label, layout) in zip(subfigure_labels,gs)
        Label(layout[1, :, Top()], label,
            fontsize = fontsize,
            padding = (-1, -1, -1, -1),
            halign = :left)
        layout.alignmode = Mixed(top=0)
    end
    """
    for g in [g_a,g_b,g_c,g_d,g_e,g_f]
        rowsize!(g, 1, Fixed(200))
    
    end
    """
    shift = 0.25

    # only put noises in the first row, unfiltered data
    f6_noise_p_scaled .= 1
    f4_noise_p_scaled .= 1
    #f3_noise_p_scaled .= 1
    
    # Define the axes and symbol data
    axes = [[ax_a, ax_b], [ax_c, ax_d], [ax_e, ax_f]]
    symbol_data = [
        [raw_harm_p_scaled, raw_noise_p_scaled, raw_entropy_scaled, artifacts],
        [f6_harm_p_scaled, f6_noise_p_scaled, f6_entropy_scaled, artifacts],
        [f4_harm_p_scaled, f4_noise_p_scaled, f4_entropy_scaled, artifacts],
        #[f3_harm_p_scaled, f3_noise_p_scaled, f3_entropy_scaled, artifacts]
    ]

    # Iterate over axes and symbol data
    for (ax12, symbols) in zip(axes, symbol_data)
        for ax in ax12
            # Scatter plot symbols for each data type
            scatter_symbol(ax, symbols[1], :rect, shift, shift)
            scatter_symbol(ax, symbols[2], :utriangle, -shift, -shift)
            scatter_symbol(ax, symbols[3], :star5, -shift, shift)

            # Plot artifact symbols
            artifact_symbol(ax, Bool.(symbols[4]), shift, -shift)

            # Add vertical and horizontal grid lines
            vlines!(ax, (0:9) .+ 0.5, color = :black)
            hlines!(ax, (0:7) .+ 0.5, color = :black)
        end
    end

    for g in [g_a,g_b,g_c,g_d,g_e,g_f]#,g_cbar,g_lab_1,g_lab_2,g_lab_3]
        colgap!(g, 0)
        rowgap!(g, 0)
        g.alignmode = Mixed(right = 0,top=0)
    end

    for g in [g_a,g_c,g_e]
        g.alignmode = Mixed(left = 0)
    end

    save(savedirname,F)

end