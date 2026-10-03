using makie_post_processing
using makie_post_processing.CairoMakie
using makie_post_processing: plot_f_unnorm_vs_vpa, input_dict_dfns,
                             check_Chodura_condition
using makie_post_processing.moment_kinetics.input_structs: boltzmann_electron_response,
                                                           boltzmann_electron_response_with_simple_sheath

function plot_f_over_vpa2(case, xmax, ymax, N_list)
    fig = Figure()
    ax = Axis(fig[1,1]; limits=(-0.5, xmax, 0.0, ymax))

    for N ∈ N_list
        if occursin("wall-plus-central-ion-source-Krook--kinetic-electrons", case)
            run_path = "$(case)-$(N)el-minimum_dt1e-7"
        else
            run_path = "$(case)-$(N)el"
        end
        ri = get_run_info(run_path; dfns=true)

        if ri.composition.electron_physics ∈ (boltzmann_electron_response,
                                              boltzmann_electron_response_with_simple_sheath)
            temp_e = nothing
        else
            temp_e = get_variable(ri, "electron_parallel_temperature")
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

        if extra_offset_upper[1,end] > 0
            vlines!(ax, cutoff_upper[1,end]; linestyle=:dash, color=l.color)
        end
    end

    Legend(fig[2,1], ax; tellwidth=false, tellheight=true)

    if occursin("kinetic-electrons", case)
        output_dir = joinpath("comparison_plots", "compare-workstation-kinetic-electrons")
        if occursin("wall", case)
            plot_name = "KE_wall_plus_central_f_over_vpa2.pdf"
        elseif occursin("central", case)
            plot_name = "KE_central_f_over_vpa2.pdf"
        elseif occursin("flat", case)
            plot_name = "KE_flat_f_over_vpa2.pdf"
        else
            error("Unrecognised kinetic electron case \"$case\"")
        end
    elseif occursin("fullf", case)
        output_dir = joinpath("comparison_plots", "compare-archer-$(basename(case))")
        plot_name = "f_over_vpa2.pdf"
    else
        output_dir = joinpath("comparison_plots", "compare-archer-$(basename(case))-mk")
        plot_name = "f_over_vpa2.pdf"
    end
    save(joinpath(output_dir, plot_name), fig)

    return fig, ax
end

function Chodura_analysis()
    plot_f_over_vpa2(joinpath("runs-archer", "central-ion-source-Krook"), 3.0, 0.6, [32, 64, 128, 256])
    plot_f_over_vpa2(joinpath("runs-archer", "flat-ion-source-Krook"), 3.0, 0.6, [32, 64, 128, 256])
    plot_f_over_vpa2(joinpath("runs-archer", "wall-plus-central-ion-source-Krook"), 3.0, 2.0, [32, 64, 128, 256])

    plot_f_over_vpa2(joinpath("runs-workstation", "central-ion-source-Krook--kinetic-electrons"), 3.0, 3.0, [32])
    plot_f_over_vpa2(joinpath("runs-workstation", "flat-ion-source-Krook--kinetic-electrons"), 3.0, 3.0, [32])
    plot_f_over_vpa2(joinpath("runs-workstation", "wall-plus-central-ion-source-Krook--kinetic-electrons"), 3.0, 50.0, [32])

    return nothing
end

Chodura_analysis()
