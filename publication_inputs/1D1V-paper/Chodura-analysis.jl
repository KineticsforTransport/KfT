using makie_post_processing
using makie_post_processing.CairoMakie
using makie_post_processing: plot_f_unnorm_vs_vpa, input_dict_dfns,
                             check_Chodura_condition
using makie_post_processing.moment_kinetics.input_structs: boltzmann_electron_response

function plot_f_over_vpa2(case, xmax, ymax)
    N_list = [32, 64, 128, 256]

    fig = Figure()
    ax = Axis(fig[1,1]; limits=(-0.5, xmax, 0.0, ymax))

    for N ∈ N_list
        ri = get_run_info(joinpath("runs-archer", "$(case)-$(N)el"); dfns=true)

        if ri.composition.electron_physics === boltzmann_electron_response
            temp_e = nothing
        else
            error("Need to get appropriate electron temperature - parallel or full "
                  * "temperature?")
        end

        f_lower = get_variable(ri, "f", iz=1)
        f_upper = get_variable(ri, "f", iz=ri.z.n_global)
        Chodura_ratio_lower, Chodura_ratio_upper, cutoff_lower, cutoff_upper,
        extra_offset_lower, extra_offset_upper =
            check_Chodura_condition(ri.r_local, ri.z_local, ri.vperp, ri.vpa,
                                    get_variable(ri, "density"),
                                    get_variable(ri, "parallel_flow"),
                                    get_variable(ri, "thermal_speed"), temp_e,
                                    ri.composition, get_variable(ri, "Er"), ri.geometry,
                                    ri.z.bc, nothing; evolve_density=ri.evolve_density,
                                    evolve_upar=ri.evolve_upar, evolve_p=ri.evolve_p,
                                    f_lower=f_lower, f_upper=f_upper,
                                    find_extra_offset=true)

        f_input = copy(input_dict_dfns["f"])
        f_input["iz0"] = ri.z.n
        l = plot_f_unnorm_vs_vpa(ri; f_over_vpa2=true, input=f_input, is=1,
                                 ax=ax, label=ri.run_name)
        #l = plot_f_unnorm_vs_vpa(ri; f_over_vpa2=true, input=f_input, is=1,
        #                         ax=ax, label=ri.run_name, scatter=true, markersize=5)

        if extra_offset_upper > 0
            vlines!(ax, cutoff_upper[1,end]; linestyle=:dash, color=l.color)
        end
    end

    Legend(fig[2,1], ax; tellwidth=false, tellheight=true)

    if occursin("fullf", case)
        output_dir = joinpath("comparison_plots", "compare-archer-$case")
    else
        output_dir = joinpath("comparison_plots", "compare-archer-$case-mk")
    end
    save(joinpath(output_dir, "f_over_vpa2.pdf"), fig)

    return fig, ax
end

function Chodura_analysis()
    plot_f_over_vpa2("central-ion-source-Krook", 3.0, 0.6)
    plot_f_over_vpa2("flat-ion-source-Krook", 3.0, 0.6)
    plot_f_over_vpa2("wall-plus-central-ion-source-Krook", 3.0, 2.0)

    return nothing
end

Chodura_analysis()
