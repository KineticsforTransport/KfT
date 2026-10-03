using makie_post_processing
using makie_post_processing: grid_points_to_faces
using makie_post_processing.CairoMakie

using Debugger
function f_2dplot_example()
    case = "wall-plus-central-ion-source-Krook"
    N = 256
    rasterize = 8
    axislabelsize = 20
    transpose_axes = true

    if occursin("fullf", case)
        output_dir = joinpath("comparison_plots", "compare-archer-$(case)")
    else
        output_dir = joinpath("comparison_plots", "compare-archer-$(case)-mk")
    end

    ri = get_run_info(joinpath("runs-archer", case * "-$(N)el"); dfns=true)

    with_theme(theme_latexfonts()) do
        fig = Figure()
        if transpose_axes
            ax = Axis(fig[1,1]; ylabel=L"v_\parallel", xlabel=L"z")
        else
            ax = Axis(fig[1,1]; xlabel=L"v_\parallel", ylabel=L"z")
        end
        ax.xlabelsize = axislabelsize
        ax.ylabelsize = axislabelsize

        plot_f_unnorm_vs_vpa_z(ri; input=nothing, electron=false, neutral=false, it=-1,
                               is=1, ax=ax, colorbar_place=fig[1,2], title="",
                               transform=identity, transpose_axes, rasterize)

        save(joinpath(output_dir, "f_2d_example.pdf"), fig)

        # Plot line showing position parallel flow - i.e. centre of w_parallel grid.
        z = ri.z.grid
        u = get_variable(ri, "parallel_flow"; it=-1, is=1, ir=1)
        if transpose_axes
            lines!(ax, z, u; color=:orange)
        else
            lines!(ax, u, z; color=:orange)
        end

        # Plot contours of logf.
        f = get_variable(ri, "f_unnorm"; it=-1, is=1, ir=1, ivperp=1)
        vpa = get_variable(ri, "vpa_unnorm"; it=-1, is=1, ir=1, ivperp=1)
        z2d = zeros(size(f))
        z2d .= reshape(z, 1, length(z))
        if transpose_axes
            contour!(ax, z2d', vpa', f';
                     levels=[10.0^i for i ∈ -16:-1],
                     color=:white, linewidth=0.5, rasterize)

        else
            contour!(ax, vpa, z2d, f;
                     levels=[10.0^i for i ∈ -16:-1],
                     color=:white, linewidth=0.5, rasterize)
        end

        save(joinpath(output_dir, "f_2d_example_fancy.pdf"), fig)
    end

    return nothing
end

f_2dplot_example()
