using makie_post_processing
using makie_post_processing.CairoMakie
using makie_post_processing: plot_f_unnorm_vs_vpa, input_dict_dfns,
                             check_Chodura_condition
using makie_post_processing.moment_kinetics.input_structs: boltzmann_electron_response,
                                                           boltzmann_electron_response_with_simple_sheath

function plot_constraints(case, N_list)
    with_theme(theme_latexfonts()) do
        fig = Figure()
        ax = Axis(fig[1,1]; yscale=log10, xlabel=L"z", ylabel="max. mean absolute coefficient")

        for N ∈ N_list
            if occursin("wall-plus-central-ion-source-Krook--kinetic-electrons", case)
                run_path = "$(case)-$(N)el-minimum_dt1e-7"
            else
                run_path = "$(case)-$(N)el"
            end
            ri = get_run_info(run_path; dfns=true)

            A = get_variable(ri, "ion_constraints_absAminus1_coefficient"; it=-1, is=1, ir=1)
            B = get_variable(ri, "ion_constraints_absB_coefficient"; it=-1, is=1, ir=1)
            C = get_variable(ri, "ion_constraints_absC_coefficient"; it=-1, is=1, ir=1)

            max_constraint = max.(A, B, C)

            lines!(ri.z.grid, max_constraint; label=L"%$N")
        end

        Legend(fig[2,1], ax; tellwidth=false, tellheight=true)

        if occursin("kinetic-electrons", case)
            output_dir = joinpath("comparison_plots", "compare-workstation-kinetic-electrons")
            if occursin("wall", case)
                plot_name = "KE_wall_plus_central_constraints.pdf"
            elseif occursin("central", case)
                plot_name = "KE_central_constraints.pdf"
            elseif occursin("flat", case)
                plot_name = "KE_flat_constraints.pdf"
            else
                error("Unrecognised kinetic electron case \"$case\"")
            end
        elseif occursin("fullf", case)
            output_dir = joinpath("comparison_plots", "compare-archer-$(basename(case))")
            if occursin("wall", case)
                plot_name = "wall_plus_central_constraints_fullf.pdf"
            elseif occursin("central", case)
                plot_name = "central_constraints_fullf.pdf"
            elseif occursin("flat", case)
                plot_name = "flat_constraints_fullf.pdf"
            else
                error("Unrecognised fullf case \"$case\"")
            end
        else
            output_dir = joinpath("comparison_plots", "compare-archer-$(basename(case))-mk")
            if occursin("wall", case)
                plot_name = "wall_plus_central_constraints.pdf"
            elseif occursin("central", case)
                plot_name = "central_constraints.pdf"
            elseif occursin("flat", case)
                plot_name = "flat_constraints.pdf"
            else
                error("Unrecognised mk case \"$case\"")
            end
        end
        save(joinpath(output_dir, plot_name), fig)
    end

    return nothing
end

function plot_constraints_all()
    plot_constraints(joinpath("runs-archer", "central-ion-source-Krook"), [32, 64, 128, 256])
    plot_constraints(joinpath("runs-archer", "flat-ion-source-Krook"), [32, 64, 128, 256])
    plot_constraints(joinpath("runs-archer", "wall-plus-central-ion-source-Krook"), [32, 64, 128, 256])

    #plot_constraints(joinpath("runs-workstation", "central-ion-source-Krook--kinetic-electrons"), [32])
    #plot_constraints(joinpath("runs-workstation", "flat-ion-source-Krook--kinetic-electrons"), [32])
    #plot_constraints(joinpath("runs-workstation", "wall-plus-central-ion-source-Krook--kinetic-electrons"), [32])

    return nothing
end

plot_constraints_all()
