using makie_post_processing
using makie_post_processing.Dates: now

function print_time_for_run(run_names...)
    run_info = get_run_info(run_names...)
    dir_name = "comparison_plots"
    open(joinpath(dir_name, "time_for_run.txt"), "a") do io
        for ri ∈ run_info
            run_time_per_step = get_variable(ri, "time_for_run_per_step")[end]
            total_time = sum(run_time_per_step)
            output_dt = get_variable(ri, "time_per_step")
            if output_dt > 0.2
                # Sim time per output step is ~1
                step_time = run_time_per_step[end]
            else
                # Sim time per output step is ~0.1
                if length(run_time_per_step) ≥ 10
                    step_time = sum(run_time_per_step[end-9:end])
                else
                    # The total run did not reach 1 unit of simulation time, so report the
                    # run time per unit simulation time for the final step, but make it
                    # negative to mark that this 'fudge' is being done.
                    step_time = - run_time_per_step[end] / step_time[end]
                end
            end
            result_string = "$(now()) $(ri.run_name) total_time=$total_time step_time=$step_time"
            println(result_string)
            println(io, result_string)
        end
    end
    return nothing
end
