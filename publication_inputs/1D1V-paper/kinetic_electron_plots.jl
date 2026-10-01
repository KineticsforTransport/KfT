using makie_post_processing
using makie_post_processing.CairoMakie

function plot_ke_case(case)
    output_directory = joinpath("comparison_plots", "compare-workstation-kinetic-electrons")

    ri = get_run_info(joinpath("runs-workstation", case); dfns=true)

    n = get_variable(ri, "density"; it=-1, is=1, ir=1)
    u = get_variable(ri, "parallel_flow"; it=-1, is=1, ir=1)
    Ti = get_variable(ri, "parallel_temperature"; it=-1, is=1, ir=1)
    ne = get_variable(ri, "electron_density"; it=-1, ir=1)
    ppar_e = get_variable(ri, "electron_parallel_pressure"; it=-1, ir=1)
    Te = ppar_e ./ ne

    if occursin("wall", case)
        prefix = "KE_wall_plus_central_"
    elseif occursin("central", case)
        prefix = "KE_central_"
    elseif occursin("flat", case)
        prefix = "KE_flat_"
    else
        error("Unrecognised case \"$case\"")
    end
    prefix = joinpath(output_directory, prefix)

    z = ri.z.grid

    fign = Figure(; size=(1200,800))
    axn = Axis(fign[1,1]; xlabel=L"z", ylabel=L"n")
    ln = lines!(axn, z, n)
    save(prefix * "density.pdf", fign)

    figu = Figure(; size=(1200,800))
    axu = Axis(figu[1,1]; xlabel=L"z", ylabel=L"u_\parallel")
    lu = lines!(axu, z, u)
    save(prefix * "parallel_flow.pdf", figu)

    figT = Figure(; size=(1200,800))
    axT = Axis(figT[1,1]; xlabel=L"z", ylabel=L"T")
    lT = lines!(axT, z, Ti; label=L"T_{i,\parallel}")
    lines!(axT, z, Te; label=L"T_{e,\parallel}")
    Legend(figT[2,1], axT; tellwidth=false, tellheight=true)
    save(prefix * "parallel_temperature.pdf", figT)

    return nothing
end

function plot_ke_all()
    for case ∈ ["central-ion-source-Krook--kinetic-electrons-32el",
                "flat-ion-source-Krook--kinetic-electrons-32el",
                "wall-plus-central-ion-source-Krook--kinetic-electrons-32el"]
        plot_ke_case(case)
    end
    return nothing
end

plot_ke_all()
