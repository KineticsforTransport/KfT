module periodic_damping

# Test generated from TOML input files

include("setup.jl")

using Base.Filesystem: tempname
using MPI

using moment_kinetics.interpolation: interpolate_to_grid_z
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info,
                                 postproc_load_variable
using moment_kinetics.utils: merge_dict_with_kwargs!

# default inputs for tests
test_input = OptionsDict(
 "output" => OptionsDict(
     "display_timing_info" => false,
     "parallel_io" => true
    ),
 "r" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1
    ),
 "vperp" => OptionsDict(
     "ngrid" => 6,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 5,
     "L" => 4.0,
     "element_spacing_option" => "coarse_tail_vperp0.5",
     "bc" => "zero"
    ),
 "evolve_moments" => OptionsDict(
     "pressure" => true,
     "density" => true,
     "moments_conservation" => true,
     "parallel_flow" => true
    ),
 "ion_species_1" => OptionsDict(
     "initial_temperature" => 2.0,
     "initial_density" => 5.0
    ),
 "krook_collisions" => OptionsDict(
     "use_krook" => true,
     "frequency_option" => "reference_parameters"
    ),
 "vpa" => OptionsDict(
     "ngrid" => 6,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 10,
     "L" => 8.0,
     "element_spacing_option" => "coarse_tails1.0",
     "bc" => "zero"
    ),
 "geometry" => OptionsDict(
     "pitch" => 1.0,
     "option" => "1D-mirror",
     "DeltaB" => 0.5
    ),
 "z" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 8,
     "L" => 13.0,
     "bc" => "periodic"
    ),
 "vpa_IC_ion_species_1" => OptionsDict(
     "initialization_option" => "gaussian",
     "density_amplitude" => 1.0,
     "temperature_amplitude" => 0.0,
     "density_phase" => 0.0,
     "upar_amplitude" => 0.0,
     "temperature_phase" => 0.0,
     "upar_phase" => 0.0
    ),
 "composition" => OptionsDict(
     "T_e" => 3.5,
     "electron_physics" => "boltzmann_electron_response",
     "n_ion_species" => 1,
     "n_neutral_species" => 0
    ),
 "ion_numerical_dissipation" => OptionsDict(
     "z_dissipation_coefficient" => 0.0
    ),
 "z_IC_ion_species_1" => OptionsDict(
     "initialization_option" => "sinusoid",
     "wavenumber" => 0,
     "density_amplitude" => 0.0,
     "temperature_amplitude" => 0.0,
     "density_phase" => 0.0,
     "upar_amplitude" => 1.0,
     "temperature_phase" => 0.0,
     "upar_phase" => 0.0
    ),
 "timestepping" => OptionsDict(
     "converged_residual_value" => 1.0e-7,
     "nstep" => 1000,
     "steady_state_residual" => true,
     "dt" => 0.05,
     "nwrite" => 500,
     "type" => "SSPRK4",
     "nwrite_dfns" => 500,
     "print_nT_live" => true
    )
)

test_input_2 = recursive_merge(test_input,
                               OptionsDict(
 "timestepping" => OptionsDict(
     "nstep" => 20000
    )
))
# Here choose the names for each test
test_input = recursive_merge(test_input,
                               OptionsDict("output" => OptionsDict("run_name" => "periodic damping Ti200 ui1")))
test_input_2 = recursive_merge(test_input_2,
                               OptionsDict("output" => OptionsDict("run_name" => "periodic damping Ti200 ui1 long")))

"""
Run a test for a single set of parameters
"""
function run_test(test_input, expected_phi; rtol=4.e-14, atol=1.e-15, args...)
    # by passing keyword arguments to run_test, args becomes a Tuple of Pairs which can be
    # used to update the default inputs

    # Make a copy to make sure nothing modifies the input Dicts defined in this test
    # script.
    input = deepcopy(test_input)

    # Convert keyword arguments to a unique name
    name = input["output"]["run_name"]
    if length(args) > 0
        name = string(name, "_", (string(k, "-", v, "_") for (k, v) in args)...)

        # Remove trailing "_"
        name = chop(name)
    end

    # Provide some progress info
    println("    - testing ", name)

    # Update default inputs with values to be changed
    merge_dict_with_kwargs!(input; args...)
    input["output"]["run_name"] = name

    # Suppress console output while running
    phi = undef
    quietoutput() do
        # run simulation
        run_moment_kinetics(input)
    end

    if global_rank[] == 0
        quietoutput() do
            # Load and analyse output
            #########################

            path = joinpath(realpath(input["output"]["base_directory"]), name)

            # open the output file(s)
            run_info = get_run_info_no_setup(path)

            # load fields data
            phi_zrt = postproc_load_variable(run_info, "phi")

            close_run_info(run_info)
            
            phi = phi_zrt[:,1,:]
        end

        # Regression test
        actual_phi = phi[begin:3:end, end]
        if expected_phi == nothing
            # Error: no expected input provided
            println("data tested would be: ", actual_phi)
            @test false
        else
            @test isapprox(actual_phi, expected_phi, rtol=rtol, atol=atol)
        end
    end
end

function runtests()
    # Create a temporary directory for test output
    test_output_directory = get_MPI_tempdir()

    @testset "periodic_damping tests" verbose=use_verbose begin
        println("periodic_damping tests")
        @testset "periodic damping Ti200 ui1" begin
            test_input["output"]["base_directory"] = test_output_directory
            run_test(test_input,
                     [5.5967425601246763, 5.5733906807391369, 5.5758764911275422, 5.5902552462771027, 5.6195252272608514, 5.6519102196310484, 5.6721369662252279, 5.6830923143176744, 5.6781836395621452, 5.6513867257101076, 5.6212679365027265])
        end
        @long @testset "periodic damping Ti200 ui1 long" begin
            test_input_2["output"]["base_directory"] = test_output_directory
            run_test(test_input_2,
                     [5.6339341024924643, 5.6339641381816268, 5.6337935534545158, 5.6332241283533602, 5.6319534019682509, 5.631111868493651, 5.6313405087209114, 5.6321743325228395, 5.633403127813235, 5.6339186669518444, 5.6339801557723295])
        end
    end
    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end

using .periodic_damping

periodic_damping.runtests()
