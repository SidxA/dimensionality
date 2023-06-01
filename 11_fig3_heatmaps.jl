#heatmap
function heatmap_numbers(savedirname)

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

#characterization map
function characterization_heatmap(savedirname)

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
        symbol = L"\xi"
        fs = 22
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

    function scatter_symbol(ax,data,symbol,x_shift,y_shift)
        symbol_sizes = [0,8,16]
        for i in 1:size(data)[1]
            str = data[i,:]
            map(enumerate(str)) do (j, c)
            scatter!(ax,i+x_shift,j+y_shift,markersize=symbol_sizes[Int(c)],marker=symbol,
            align = (:center,:center),
            color = :black,
            )

            end
        end
    end



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


    Colorbar(g_cbar[1,1:6], hm_a, vertical = false)
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

    shift = 0.25

    axes = [[ax_a,ax_b],[ax_c,ax_d],[ax_e,ax_f],[ax_g,ax_h]]
    symbol_data = [[raw_harm_p_scaled,raw_noise_p_scaled,raw_entropy_scaled,artifacts],
    [f6_harm_p_scaled,f6_noise_p_scaled,f6_entropy_scaled,artifacts],
    [f4_harm_p_scaled,f4_noise_p_scaled,f4_entropy_scaled,artifacts],
    [f3_harm_p_scaled,f3_noise_p_scaled,f3_entropy_scaled,artifacts]]

    for (ax12,symbols) in zip(axes,symbol_data)
        for ax in ax12
            scatter_symbol(ax,symbols[1],:rect,shift,shift)
            scatter_symbol(ax,symbols[2],:utriangle,-shift,-shift)
            scatter_symbol(ax,symbols[3],:star6,-shift,shift)
            artifact_symbol(ax,Bool.(symbols[4]),shift,-shift)

            vlines!(ax,(0:9) .+ 0.5,color = :black)
            hlines!(ax,(0:6) .+ 0.5,color = :black)
        end

    end

    save(savedirname,F)

end