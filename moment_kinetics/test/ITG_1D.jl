"""
Module for testing the 1D ITG helical field case.

After running the simulation, the growth rate is calculated for the 
single initialised mode using a DFFT or something.

Growth rate in code units (check page 26 of FP4):

γ = k_z * ρ_i * v_ti/(2*sqrt(L_B * L_T)) - D * k_z^2

in units in the code, this is 

γ = k_z/(sqrt(L_B * L_T)) - D * k_z^2

"""

module ITG_1D

# Test generated from TOML input files

include("setup.jl")

using Base.Filesystem: tempname
using MPI

using moment_kinetics.interpolation: interpolate_to_grid_z
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info,
                                 postproc_load_variable, get_variable
using moment_kinetics.utils: merge_dict_with_kwargs!
using moment_kinetics.analysis: get_growth_rate_of_box_mode_1D

# default inputs for tests
test_input = OptionsDict(
 "output" => OptionsDict(
     "display_timing_info" => false
    ),
 "r" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1
    ),
 "vperp" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 5,
     "L" => 6.0,
     "bc" => "zero"
    ),
 "evolve_moments" => OptionsDict(
     "pressure" => false,
     "density" => false,
     "parallel_flow" => false
    ),
 "ion_species_1" => OptionsDict(
     "initial_temperature" => 1.0,
     "initial_density" => 1.0,
     "L_T" => -0.005
    ),
 "krook_collisions" => OptionsDict(
     "use_krook" => true,
     "frequency_option" => "reference_parameters"
    ),
 "vpa" => OptionsDict(
     "ngrid" => 6,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 10,
     "L" => 9.0,
     "bc" => "zero"
    ),
 "geometry" => OptionsDict(
     "pitch" => 0.0,
     "option" => "1D-Helical-ITG",
     "rhostar" => 1.0,
     "dBdr_constant" => 0.05
    ),
 "z" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 16,
     "L" => 1.0,
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
     "T_e" => 1.0,
     "electron_physics" => "boltzmann_electron_response",
     "n_ion_species" => 1,
     "n_neutral_species" => 0
    ),
 "ion_numerical_dissipation" => OptionsDict(
     "z_dissipation_coefficient" => 0.1,
     "z_dissipation_degree" => 2
    ),
 "z_IC_ion_species_1" => OptionsDict(
     "initialization_option" => "sinusoid",
     "wavenumber" => 1.0,
     "density_amplitude" => 1.0e-6,
     "temperature_amplitude" => 0.0,
     "density_phase" => 0.0,
     "upar_amplitude" => 0.0,
     "temperature_phase" => 0.0,
     "upar_phase" => 0.0
    ),
 "timestepping" => OptionsDict(
     "nstep" => 960,
     "steady_state_residual" => true,
     "dt" => 0.0005,
     "nwrite" => 20,
     "type" => "SSPRK4",
     "nwrite_dfns" => 20,
     "print_nT_live" => true
    )
)

test_input_2 = recursive_merge(test_input,
                               OptionsDict(
 "ion_species_1" => OptionsDict(
     "L_T" => -0.01
    )
))
test_input_3 = recursive_merge(test_input,
                               OptionsDict(
 "ion_species_1" => OptionsDict(
     "L_T" => -0.015
    )
))
# Here choose the names for each test
test_input = recursive_merge(test_input,
                               OptionsDict("output" => OptionsDict("run_name" => "L_B 20, L_T 0.005")))
test_input_2 = recursive_merge(test_input_2,
                               OptionsDict("output" => OptionsDict("run_name" => "L_B 20, L_T 0.01")))
test_input_3 = recursive_merge(test_input_3,
                               OptionsDict("output" => OptionsDict("run_name" => "L_B 20, L_T 0.015")))

"""
Run a test for a single set of parameters
"""
function run_test(test_input, expected_phi, growth_rate; rtol=5.e-6, atol=1.e-15, args...)
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
    phi = nothing
    γ = nothing
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
            phi = phi_zrt[:,1,:]

            γ = get_growth_rate_of_box_mode_1D(run_info)

            close_run_info(run_info)
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

        # test that growth rate of box scale mode was within tolerance of original
        @test isapprox(γ, growth_rate, rtol=rtol, atol=atol)

    end
end

function runtests()
    # Create a temporary directory for test output
    test_output_directory = get_MPI_tempdir()

    @testset "ITG_1D tests" verbose=use_verbose begin
        println("ITG_1D tests")
        @testset "try_large_LB_1" begin
            test_input["output"]["base_directory"] = test_output_directory
            run_test(test_input,
                     [-0.0009912671086371701, -0.0010578445129206591,
                      -0.0010277748655096322, -0.0009261378528561237,
                      -0.0006909340108779558, -0.0003672088399582509,
                      -8.929205139622209e-5, 0.00019253598120273565, 0.000533844122121934,
                      0.0008016700310614562, 0.0009378552700426428, 0.0010095313054974608,
                      0.0010042731993665632, 0.0008935672276749854, 0.0007495910898176804,
                      0.0005640264268444433, 0.000273200894806827, -4.55948749131495e-5,
                      -0.0002885949913970617, -0.0005188627903846614,
                      -0.0007811744731474714, -0.0009686449462972152],
                     [15.910306768328583])
        end
        @testset "try_large_LB_2" begin
            test_input_2["output"]["base_directory"] = test_output_directory
            run_test(test_input_2,
                     [-7.485568887933946e-5, -7.718810515670286e-5, -7.425227088535807e-5,
                      -6.74420629642205e-5, -5.3279525623125804e-5, -3.471649855738671e-5,
                      -1.898580575720367e-5, -2.929741594114896e-6, 1.7018148812146337e-5,
                      3.350203389897328e-5, 4.253663178141189e-5, 4.794167222826572e-5,
                      4.8974137351331544e-5, 4.3010919073798413e-5, 3.4225254828548317e-5,
                      2.2365397683751354e-5, 3.3047176783339036e-6,
                      -1.7663983823560278e-5, -3.3396340594388e-5, -4.7888988421652556e-5,
                      -6.357459071200297e-5, -7.375571471966637e-5],
                     [10.07945303817917])
        end
        @testset "try_large_LB_3" begin
            test_input_3["output"]["base_directory"] = test_output_directory
            run_test(test_input_3,
                     [-3.1630698110819085e-5, -3.2262694954050406e-5,
                      -3.137955278493604e-5, -2.9378264634412e-5, -2.5248879242288195e-5,
                      -1.9859763129996062e-5, -1.530365183161917e-5,
                      -1.0660075072958585e-5, -4.89856276473545e-6,
                      -1.4451754438250728e-7, 2.455855931919879e-6, 4.0046117011427155e-6,
                      4.280679021136262e-6, 2.5287779315793404e-6, -3.5820955177484904e-8,
                      -3.4910757157391913e-6, -9.036115759536784e-6,
                      -1.5126462980064575e-5, -1.96879147593425e-5, -2.388139248050448e-5,
                      -2.8404556903081285e-5, -3.131905945445779e-5],
                     [7.487059462842807])
        end
    end
    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end

using .ITG_1D

ITG_1D.runtests()
