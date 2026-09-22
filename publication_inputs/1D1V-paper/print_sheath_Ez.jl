using makie_post_processing
using makie_post_processing.Dates: now
using makie_post_processing.StatsBase: mean

function print_sheath_Ez(run_names...)
    run_info = get_run_info(run_names...)
    dir_name = "comparison_plots"
    open(joinpath(dir_name, "sheath_Ez.txt"), "a") do io
        for ri ∈ run_info
            Ez = get_variable(ri, "Ez"; it=-1, ir=1)
            sheath_Ez = mean(abs.(Ez[[1,end]]))
            result_string = "$(now()) $(ri.run_name) $sheath_Ez"
            println(result_string)
            println(io, result_string)
        end
    end
    return nothing
end
