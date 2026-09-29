using Dates
using LaTeXStrings
using makie_post_processing
using makie_post_processing: regrid_variable, _run_info_to_coords
using makie_post_processing.CairoMakie

function get_error(ri, ri_highres, variable_name::String, atol::Number)
    v_highres = get_variable(ri_highres, variable_name; it=-1)
    v_before_interp = get_variable(ri, variable_name; it=-1)
    v = regrid_variable(v_before_interp, _run_info_to_coords(ri_highres),
                        _run_info_to_coords(ri), ri_highres, ri.evolve_density,
                        ri.evolve_upar, ri.evolve_p)

    # Calculate a relative error, with absolute tolerance.
    err = @. (v_highres - v) / (0.5 * (abs(v_highres) + abs(v)) + atol)
    println("      max differences ", extrema(v_highres .- v))
    minind = argmin(err)

    error_rms = sqrt(sum(err.^2) / length(err))
    error_max = maximum(abs.(err))

    #return error_rms, error_max
    return error_rms, error_max, err
end

function plot_error_scaling(case::String, variable_name::String, ax_rms, ax_max,
                            highres_N::Number, lowres_N_list::AbstractVector,
                            output_dir::String; atol=1.0e-12)
    println("  $variable_name")

    function get_name(N)
        return joinpath("runs-archer", case * "-$(N)el")
    end

    ri_highres = redirect_stdout(devnull) do
        get_run_info(get_name(highres_N); dfns=true)
    end

    error_rms_list = Float64[]
    error_max_list = Float64[]
    open(joinpath(output_dir, "relative_errors.txt"), "a") do io
        for N ∈ lowres_N_list
            ri =  redirect_stdout(devnull) do
                get_run_info(get_name(N); dfns=true)
            end
            #error_rms, error_max = get_error(ri, ri_highres, variable_name, atol)
            error_rms, error_max, err = get_error(ri, ri_highres, variable_name, atol)
            println("    N=$N error_rms=$error_rms error_max=$error_max")
            push!(error_rms_list, error_rms)
            push!(error_max_list, error_max)

            println(io, "$(now()) $(ri.run_name) $variable_name atol=$atol N=$N error_rms=$error_rms error_max=$error_max")
        end
    end

    # Plots for rms errors.
    #######################

    # Plot simulation errors
    s = scatter!(ax_rms, lowres_N_list, error_rms_list; marker=:x, label=variable_name)

    # Plot N^(-2) scaling to see if it fits.
    scaling_line_rms = @. error_rms_list[1] * (lowres_N_list[1] / lowres_N_list)^2
    lines!(ax_rms, lowres_N_list, scaling_line_rms; color=s.color, linestyle=:dot)

    # Plots for max errors.
    #######################

    # Plot simulation errors
    s = scatter!(ax_max, lowres_N_list, error_max_list; marker=:x, label=variable_name)

    # Plot N^(-2) scaling to see if it fits.
    scaling_line_max = @. error_max_list[1] * (lowres_N_list[1] / lowres_N_list)^2
    lines!(ax_max, lowres_N_list, scaling_line_max; color=s.color, linestyle=:dot)

    return nothing
end

function make_convergence_plots(case)
    println(case)
    highres_N = 256
    lowres_N_list = [32, 64, 128]
    if occursin("fullf", case)
        output_dir = joinpath("comparison_plots", "compare-archer-$case")
    else
        output_dir = joinpath("comparison_plots", "compare-archer-$(case)-mk")
    end

    with_theme(theme_latexfonts()) do
        fig_rms = Figure()
        ax_rms = Axis(fig_rms[1,1], yscale=log10, xscale=log2, xticks=2.0.^lowres_N_list, xlabel=L"N", ylabel="RMS relative error")
        fig_max = Figure()
        ax_max = Axis(fig_max[1,1], yscale=log10, xscale=log2, xticks=2.0.^lowres_N_list, xlabel=L"N", ylabel="Maximum absolute error")

        #for variable_name ∈ ["density", "parallel_flow", "parallel_temperature", "f"]
        for variable_name ∈ ["density", "parallel_flow", "parallel_temperature"]
            plot_error_scaling(case, variable_name, ax_rms, ax_max, highres_N,
                               lowres_N_list, output_dir)
        end

        Legend(fig_rms[2,1], ax_rms; tellwidth=false, tellheight=true)
        save(joinpath(output_dir, "error_rms.pdf"), fig_rms)

        Legend(fig_max[2,1], ax_max; tellwidth=false, tellheight=true)
        save(joinpath(output_dir, "error_max.pdf"), fig_max)
    end
end

function make_all_convergence_plots()
    make_convergence_plots("wall-plus-central-ion-source-Krook")
    return nothing
end

make_all_convergence_plots()
