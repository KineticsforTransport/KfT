module Multi_Source_Tests

# Test for multi source functionality 

include("setup.jl")

using Base.Filesystem: tempname
using MPI

using moment_kinetics.interpolation: interpolate_to_grid_z
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info,
                                 postproc_load_variable
using moment_kinetics.utils: merge_dict_with_kwargs!

# default inputs for tests
test_input = OptionsDict("composition" => OptionsDict("n_ion_species" => 1,
                                                      "n_neutral_species" => 0,
                                                      "electron_physics" => "boltzmann_electron_response",
                                                      "T_e" => 0.8),
                         "ion_species_1" => OptionsDict("initial_density" => 5.0,
                                                        "initial_temperature" => 0.3333333333333333),
                         "z_IC_ion_species_1" => OptionsDict("initialization_option" => "gaussian",
                                                             "density_amplitude" => 0.001,
                                                             "density_phase" => 0.0,
                                                             "upar_amplitude" => 1.4142135623730951,
                                                             "upar_phase" => 0.0,
                                                             "temperature_amplitude" => 0.0,
                                                             "temperature_phase" => 0.0),
                         "vpa_IC_ion_species_1" => OptionsDict("initialization_option" => "gaussian",
                                                               "density_amplitude" => 1.0,
                                                               "density_phase" => 0.0,
                                                               "upar_amplitude" => 0.0,
                                                               "upar_phase" => 0.0,
                                                               "temperature_amplitude" => 0.0,
                                                               "temperature_phase" => 0.0),
                         "output" => OptionsDict("run_name" => "Multi Source Test",
                                                 "display_timing_info" => false),
                         "evolve_moments" => OptionsDict("density" => true,
                                                         "parallel_flow" => true,
                                                         "pressure" => true),
                         "timestepping" => OptionsDict("nstep" => 200,
                                                       "dt" => 0.0007071067811865475,
                                                       "nwrite" => 100,
                                                       "type" => "SSPRK4",
                                                       "nwrite_dfns" => 100,
                                                       "steady_state_residual" => true),
                         "r" => OptionsDict("ngrid" => 1,
                                            "nelement" => 1),
                         "z" => OptionsDict("ngrid" => 5,
                                            "nelement" => 20,
                                            "bc" => "wall",
                                            "discretization" => "chebyshev_pseudospectral"),
                         "vpa" => OptionsDict("ngrid" => 6,
                                              "nelement" => 20,
                                              "L" => 27.712812921102035,
                                              "bc" => "zero",
                                              "discretization" => "chebyshev_pseudospectral",
                                              "element_spacing_option" => "coarse_tails8.660254037844386"),
                         "ion_source_1" => OptionsDict("active" => true,
                                                       "z_profile" => "super_gaussian_4",
                                                       "z_width" => 0.275816,
                                                       "source_strength" => 7.0710678118654755,
                                                       "source_T" => 1.5,
                                                       "source_type" => "Maxwellian"),
                         "ion_source_2" => OptionsDict("active" => true,
                                                       "z_profile" => "wall_exp_decay",
                                                       "z_width" => 0.15,
                                                       "source_strength" => 7.0710678118654755,
                                                       "source_T" => 0.1,
                                                       "source_type" => "Maxwellian"),
                         "krook_collisions" => OptionsDict("use_krook" => true,
                                                           "frequency_option" => "reference_parameters"),
                        )

if global_size[] > 2 && global_size[] % 2 == 0
    # Test using distributed-memory
    test_input["z"]["nelement_local"] = test_input["z"]["nelement"] ÷ 2
end

test_input_1 = recursive_merge(test_input,
                               OptionsDict("output" => OptionsDict("run_name" => "two_ion_sources_moments")))
test_input_2 = recursive_merge(test_input_1,
                               OptionsDict("output" => OptionsDict("run_name" => "two_ion_sources"),
                                           "evolve_moments" => OptionsDict("density" => false,
                                                                           "parallel_flow" => false,
                                                                           "pressure" => false),
                                           "vpa" => OptionsDict("L" => 22.627416997969522,
                                                                "element_spacing_option" => "coarse_tails7.0710678118654755"),
                                          ))
test_input_3 = recursive_merge(test_input,
                               OptionsDict("output" => OptionsDict("run_name" => "PI_controller_sources_moments"),
                                           "ion_source_1" => OptionsDict("active" => true,
                                                                         "z_profile" => "super_gaussian_4",
                                                                         "z_width" => 0.275816,
                                                                         "source_strength" => 7.0710678118654755,
                                                                         "source_T" => 1.5,
                                                                         "source_type" => "temperature_midpoint_control",
                                                                         "PI_temperature_controller_I" => 30.0,
                                                                         "PI_temperature_controller_P" => 21.213203435596427,
                                                                         "PI_temperature_target_amplitude" => 0.3333333333333333),
                                           "ion_source_2" => OptionsDict("active" => true,
                                                                         "z_profile" => "wall_exp_decay",
                                                                         "z_width" => 0.15,
                                                                         "source_strength" => 7.0710678118654755,
                                                                         "source_T" => 0.1,
                                                                         "source_type" => "density_midpoint_control",
                                                                         "PI_density_controller_I" => 2.0,
                                                                         "PI_density_controller_P" => 1.4142135623730951,
                                                                         "PI_density_target_amplitude" => 1.15)))
test_input_4 = recursive_merge(test_input,
                               OptionsDict("output" => OptionsDict("run_name" => "PI_controller_sources"),
                                           "ion_source_1" => OptionsDict("active" => true,
                                                                         "z_profile" => "super_gaussian_4",
                                                                         "z_width" => 0.275816,
                                                                         "source_strength" => 7.0710678118654755,
                                                                         "source_T" => 1.5,
                                                                         "source_type" => "temperature_midpoint_control",
                                                                         "PI_temperature_controller_I" => 30.0,
                                                                         "PI_temperature_controller_P" => 21.213203435596427,
                                                                         "PI_temperature_target_amplitude" => 0.3333333333333333),
                                           "ion_source_2" => OptionsDict("active" => true,
                                                                         "z_profile" => "wall_exp_decay",
                                                                         "z_width" => 0.15,
                                                                         "source_strength" => 7.0710678118654755,
                                                                         "source_T" => 0.1,
                                                                         "source_type" => "density_midpoint_control",
                                                                         "PI_density_controller_I" => 2.0,
                                                                         "PI_density_controller_P" => 1.4142135623730951,
                                                                         "PI_density_target_amplitude" => 1.15),
                                           "evolve_moments" => OptionsDict("density" => false,
                                                                           "parallel_flow" => false,
                                                                           "pressure" => false),
                                           "vpa" => OptionsDict("L" => 22.627416997969522,
                                                                "element_spacing_option" => "coarse_tails7.0710678118654755"),
                                          ))
test_input_3["timestepping"] = recursive_merge(test_input_3["timestepping"],
                                               OptionsDict("nstep" => 500))
test_input_4["timestepping"] = recursive_merge(test_input_4["timestepping"],
                                               OptionsDict("nstep" => 500))


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

    @testset "Multi Source Tests" verbose=use_verbose begin
        println("multi source tests")
        @testset "multi source test 1" begin
            test_input_1["output"]["base_directory"] = test_output_directory
            run_test(test_input_1,
                     [1.1747328170890627, 1.246347218686613, 1.2305695288471379,
                      1.2028168291962527, 1.1702014355384738, 1.1471225228492226,
                      1.1381555217443464, 1.133938810259904, 1.1329172796177474,
                      1.13488256354748, 1.134611494149954, 1.1332181647750295,
                      1.1320162768837105, 1.1311623022422361, 1.1313127824602947,
                      1.1321662870941431, 1.133542411012697, 1.1349936375290441,
                      1.1341722835743837, 1.1327719236839753, 1.1346110710325858,
                      1.142234725869274, 1.1550677772969653, 1.17524334540261,
                      1.2093593042554551, 1.2433925046880663, 1.2399484673075234])
        end
        @testset "multi source test 2" begin
            test_input_2["output"]["base_directory"] = test_output_directory
            run_test(test_input_2,
                     [1.1655594865309673, 1.2495047470884797, 1.2309976172712136,
                      1.2034822179184976, 1.1697999553188756, 1.1463454370580293,
                      1.136717345272336, 1.1331639419699535, 1.1336775694711045,
                      1.1351193814970546, 1.1350547486105123, 1.1339927128368603,
                      1.1320762538074045, 1.131007894606059, 1.1312740709930746,
                      1.1323814549394349, 1.1342953190135052, 1.1352369362109496,
                      1.1346639276500559, 1.1333979576182693, 1.1335446883638753,
                      1.141136046855272, 1.1545894547758782, 1.174956401346126,
                      1.2097977104198447, 1.242751958455989, 1.2408823667463023])
        end
        @testset "multi source test 3" begin
            test_input_3["output"]["base_directory"] = test_output_directory
            run_test(test_input_3,
                     [0.6854986995190636, 0.7943798644013363, 0.7861596518067642,
                      0.7505238849965936, 0.681698346478118, 0.599231996914371,
                      0.535298000658603, 0.4743742423283961, 0.4003574060130213,
                      0.3352590510030524, 0.2918150980898613, 0.2538893921397605,
                      0.2170516870657307, 0.2013722716247333, 0.2049935346615851,
                      0.2219806203170459, 0.26190931672074985, 0.31495102632956107,
                      0.361148208188496, 0.4124201697723723, 0.4878733440007258,
                      0.5699909420226217, 0.6343362820709264, 0.6947443541178497,
                      0.7602641737195206, 0.7962315771583066, 0.7757122323776903])
        end
        @testset "multi source test 4" begin
            test_input_4["output"]["base_directory"] = test_output_directory
            run_test(test_input_4,
                     [0.6443315132621561, 0.7990871379202669, 0.7872536361051092,
                      0.7511513953809122, 0.6820933461434059, 0.5997816863381101,
                      0.5356329803704418, 0.47410191223764553, 0.39955191229989895,
                      0.3348497305247017, 0.2915175401928749, 0.2540607682538532,
                      0.2171685031107371, 0.20023433504118754, 0.20439968324649405,
                      0.2222271521839152, 0.26196073466944214, 0.31457305775480976,
                      0.3605894281797814, 0.41161026298505166, 0.487748641518033,
                      0.5705349226568462, 0.6348381347415498, 0.6951380823044662,
                      0.7608863511182651, 0.7972185332770916, 0.7792981485288654])
        end
    end
    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end

using .Multi_Source_Tests

Multi_Source_Tests.runtests()
