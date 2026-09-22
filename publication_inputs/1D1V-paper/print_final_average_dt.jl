using makie_post_processing
using makie_post_processing.Dates: now

function print_final_average_dt(run_names...)
    run_info = get_run_info(run_names...)
    dir_name = "comparison_plots"
    open(joinpath(dir_name, "final_average_dt.txt"), "a") do io
        for ri ∈ run_info
            dt = get_variable(ri, "average_successful_dt")[end]
            result_string = "$(now()) $(ri.run_name) $dt"
            println(result_string)
            println(io, result_string)
        end
    end
    return nothing
end
