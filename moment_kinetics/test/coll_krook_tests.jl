module coll_krook_tests

# Test generated from TOML input files

include("setup.jl")

using Base.Filesystem: tempname
using MPI

using moment_kinetics.interpolation: interpolate_to_grid_z
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info,
                                 postproc_load_variable
using moment_kinetics.utils: merge_dict_with_kwargs!

# default inputs for tests
dkions_n10 = OptionsDict(
 "output" => OptionsDict(
     "display_timing_info" => false
    ),
 "r" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1
    ),
 "evolve_moments" => OptionsDict(
     "pressure" => true,
     "density" => true,
     "moments_conservation" => true,
     "parallel_flow" => true
    ),
 "ion_species_1" => OptionsDict(
     "initial_temperature" => 0.3333333333333333,
     "initial_density" => 40.0
    ),
 "krook_collisions" => OptionsDict(
     "use_krook" => true,
     "frequency_option" => "reference_parameters"
    ),
 "vpa" => OptionsDict(
     "ngrid" => 6,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 31,
     "L" => 14.0,
     "element_spacing_option" => "coarse_tails8.660254037844386",
     "bc" => "zero"
    ),
 "z" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "gausslegendre_pseudospectral",
     "nelement" => 20,
     "L" => 1.4,
     "bc" => "wall"
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
     "T_e" => 0.8,
     "electron_physics" => "boltzmann_electron_response",
     "n_ion_species" => 1,
     "n_neutral_species" => 0
    ),
 "ion_source_2" => OptionsDict(
     "source_type" => "density_midpoint_control",
     "source_T" => 0.1,
     "active" => true,
     "PI_density_controller_I" => 5.0,
     "source_strength" => 40.0,
     "z_profile" => "wall_exp_decay",
     "PI_density_controller_P" => 7.0,
     "PI_density_target_amplitude" => 10.0,
     "z_width" => 0.2
    ),
 "z_IC_ion_species_1" => OptionsDict(
     "initialization_option" => "gaussian",
     "density_amplitude" => 0.001,
     "temperature_amplitude" => 0.0,
     "density_phase" => 0.0,
     "upar_amplitude" => 1.4142135623730951,
     "temperature_phase" => 0.0,
     "upar_phase" => 0.0
    ),
 "ion_source_1" => OptionsDict(
     "PI_temperature_controller_I" => 500.0,
     "source_type" => "temperature_midpoint_control",
     "source_T" => 1.5,
     "active" => true,
     "PI_temperature_target_amplitude" => 0.3333333333333333,
     "source_strength" => 14.0,
     "z_profile" => "super_gaussian_4",
     "PI_temperature_controller_P" => 500.0,
     "z_width" => 0.38
    ),
 "timestepping" => OptionsDict(
     "nstep" => 4000,
     "steady_state_residual" => true,
     "dt" => 0.0005,
     "nwrite" => 2000,
     "type" => "SSPRK4",
     "nwrite_dfns" => 2000,
     "print_nT_live" => true
    )
)

coll_krook_n10 = recursive_merge(dkions_n10,
                               OptionsDict(
 "vpa" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1,
    ),
 "composition" => OptionsDict(
     "ion_physics" => "coll_krook_ions"
    ),
 "z" => OptionsDict(
     "ngrid" => 6,
     "nelement" => 40,
    ),
 "timestepping" => OptionsDict(
     "nstep" => 25000,
     "dt" => 2.0e-5,
     "nwrite" => 2500,
     "nwrite_dfns" => 2500,
    )
))

dkions_n1 = recursive_merge(dkions_n10,
                               OptionsDict(
 "ion_source_2" => OptionsDict(
     "source_strength" => 4.0,
     "PI_density_target_amplitude" => 1.0
    ),
 "vpa" => OptionsDict(
     "nelement" => 41
    ),
 "ion_source_1" => OptionsDict(
     "source_strength" => 1.4
    ),
 "timestepping" => OptionsDict(
     "nstep" => 400,
     "nwrite" => 200,
     "nwrite_dfns" => 200
    ),
 "z" => OptionsDict(
     "nelement" => 32
    ),
 "ion_species_1" => OptionsDict(
     "initial_density" => 4.0
    )
))

coll_krook_n1 = recursive_merge(dkions_n1,
                               OptionsDict(
 "vpa" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1,
    ),
 "z" => OptionsDict(
     "ngrid" => 6,
     "nelement" => 40,
    ),
 "composition" => OptionsDict(
     "ion_physics" => "coll_krook_ions"
    ),
 "timestepping" => OptionsDict(
     "nstep" => 250000,
     "dt" => 2.0e-6,
     "nwrite" => 25000,
     "nwrite_dfns" => 25000,
    )
))

dkions_n100 = recursive_merge(dkions_n10,
                               OptionsDict(
 "ion_source_2" => OptionsDict(
     "PI_density_controller_I" => 500.0,
     "source_strength" => 200.0,
     "PI_density_controller_P" => 700.0,
     "PI_density_target_amplitude" => 100.0
    ),
 "vpa" => OptionsDict(
     "nelement" => 41
    ),
 "ion_source_1" => OptionsDict(
     "source_strength" => 140.0
    ),
 "timestepping" => OptionsDict(
     "nstep" => 80,
     "dt" => 0.0002,
     "nwrite" => 40,
     "nwrite_dfns" => 40
    ),
 "z" => OptionsDict(
     "nelement" => 32
    ),
 "ion_species_1" => OptionsDict(
     "initial_density" => 130.0
    )
))

coll_krook_n100 = recursive_merge(dkions_n100,
                               OptionsDict(
 "vpa" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1,
    ),
 "z" => OptionsDict(
     "ngrid" => 6,
     "nelement" => 40,
    ),
 "composition" => OptionsDict(
     "ion_physics" => "coll_krook_ions"
    ),
 "timestepping" => OptionsDict(
     "nstep" => 25000,
     "dt" => 2.0e-5,
     "nwrite" => 2500,
     "nwrite_dfns" => 2500,
    )
))


# Here choose the names for each test
dkions_n10 = recursive_merge(dkions_n10,
                               OptionsDict("output" => OptionsDict("run_name" => "dkions_n10.0_for_test_generation")))
coll_krook_n10 = recursive_merge(coll_krook_n10,
                               OptionsDict("output" => OptionsDict("run_name" => "coll_krook_n10.0")))
dkions_n1 = recursive_merge(dkions_n1,
                               OptionsDict("output" => OptionsDict("run_name" => "dkions_n1.0_for_test_generation")))
coll_krook_n1 = recursive_merge(coll_krook_n1,
                               OptionsDict("output" => OptionsDict("run_name" => "coll_krook_n1.0")))
dkions_n100 = recursive_merge(dkions_n100,
                               OptionsDict("output" => OptionsDict("run_name" => "dkions_n100.0_for_test_generation")))
coll_krook_n100 = recursive_merge(coll_krook_n100,
                               OptionsDict("output" => OptionsDict("run_name" => "coll_krook_n100.0")))

"""
Run a test for a single set of parameters
"""
function run_test_with_restart(dkions_test_input, coll_krook_test_input, expected_phi; rtol=1.e-12, atol=1.e-12, args...)
    # by passing keyword arguments to run_test, args becomes a Tuple of Pairs which can be
    # used to update the default inputs

    # Make a copy to make sure nothing modifies the input Dicts defined in this test
    # script.
    dkions_input = deepcopy(dkions_test_input)
    coll_krook_input = deepcopy(coll_krook_test_input)

     # Convert keyword arguments to a unique name
    dkions_name = dkions_input["output"]["run_name"]
    if length(args) > 0
        dkions_name = string(dkions_name, "_", (string(k, "-", v, "_") for (k, v) in args)...)  
        # Remove trailing "_"
        dkions_name = chop(dkions_name)
    end

    # Convert keyword arguments to a unique name
    coll_krook_name = coll_krook_input["output"]["run_name"]
    if length(args) > 0
        coll_krook_name = string(coll_krook_name, "_", (string(k, "-", v, "_") for (k, v) in args)...)

        # Remove trailing "_"
        coll_krook_name = chop(coll_krook_name)
    end

    # Provide some progress info
    println("    - testing ", coll_krook_name)

    # Update default inputs with values to be changed
    merge_dict_with_kwargs!(dkions_input; args...)
    dkions_input["output"]["run_name"] = dkions_name
    # Suppress console output while running
    phi = undef
    quietoutput() do
        # run simulation
        run_moment_kinetics(dkions_input)

        # now run the coll_krook simulation restarting from the dkions output
        name_of_restart_file = dkions_name * ".dfns.h5"
        run_moment_kinetics(coll_krook_input, restart = joinpath(
            realpath(dkions_input["output"]["base_directory"]),
            dkions_name, name_of_restart_file))
    end

    if global_rank[] == 0
        quietoutput() do
            # Load and analyse output
            #########################

            path = joinpath(realpath(coll_krook_input["output"]["base_directory"]), coll_krook_name)

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
    test_output_directory = "temptest" #get_MPI_tempdir()

    @testset "coll_krook tests" verbose=use_verbose begin
        println("coll_krook tests")
        @testset "coll_krook_test_n1.0" begin
            dkions_n1["output"]["base_directory"] = test_output_directory
            coll_krook_n1["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n1, coll_krook_n1,
                        [-0.18990806457834905, -0.040411399053795705, -0.0059686171914708425,
                       0.02647678521153092, 0.03353347837077063, 0.02843686230638412,
                       0.02320856111212149, 0.018120698656462857, 0.009313431684651187,
                       0.003634592793196831, -0.0037721362226705263, -0.010588983343425075,
                       -0.014934667423757698, -0.020340874395825895, -0.02257199183477208,
                       -0.02414658081472189, -0.0240250727053883, -0.022925847193656396,
                       -0.019611860457376368, -0.01677043450209867, -0.012215390197868882,
                       -0.00711516970114661, -0.0031598054032586026, 0.00325690779074695,
                       0.007145374484091978, 0.012150817370578517, 0.016760525763564817,
                       0.01984908299610518, 0.024171424249091716, 0.026419321646314383,
                       0.028903682545368528, 0.030724917514059848, 0.0316283075834402,
                       0.03227637920894838, 0.03218027999445438, 0.031440050280496186,
                       0.030000213185584875, 0.02849780555867322, 0.02532683154295361,
                       0.02290633001414765, 0.019113178609163277, 0.014769978673634368,
                       0.011263053654716897, 0.005196136561357638, 0.0012551095805596983,
                       -0.004146259712395947, -0.009428658361320522, -0.013098264083817895,
                       -0.018271294869465878, -0.02083517031576956, -0.0232705160848192,
                       -0.024281352448120547, -0.023982751563879908, -0.02155523583577422,
                       -0.018873137952978125, -0.013934033603235103, -0.007664364828634882,
                       -0.0024784783864515263, 0.006500829755695582, 0.01214091882828509,
                       0.01941143792990163, 0.026044582712643746, 0.029373654960129877,
                       0.032244659248226684, 0.018310346973676683, -0.013465210889757849,
                       -0.07242589162733908])
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                            [1.5839307267261333, 1.8460266737892335, 1.9053605594573282,
                            1.9460747875423872, 1.9554226822621759, 1.950374100328963,
                            1.9360223413695383, 1.920945386430461, 1.8934464543089184,
                            1.8750317487332673, 1.8507452581875456, 1.8270054158099018,
                            1.8105294711050708, 1.785734358682708, 1.7717646620464595,
                            1.7542733360382567, 1.7387738504071526, 1.728340641805038,
                            1.7135722902588437, 1.7055122036320765, 1.6961576676138421,
                            1.6883033237069478, 1.6836158927692424, 1.6772954408439678,
                            1.673905891505542, 1.6694266806562692, 1.6654587304533828,
                            1.662735080287101, 1.6590097852090988, 1.6570242702448006,
                            1.6549313888286001, 1.6532972676615498, 1.6525560053283705,
                            1.6519578852790222, 1.6520747347664906, 1.6526519259272885,
                            1.6539548957716044, 1.655215634233552, 1.6580072772193337,
                            1.6600747184523905, 1.663440249207427, 1.667169103137065,
                            1.670273439075196, 1.675583727409746, 1.6791655991980052,
                            1.6846391584335305, 1.6916038987324962, 1.6976766428093133,
                            1.70949384466913, 1.7179205997477494, 1.73090633241958,
                            1.7454282981972344, 1.7574173312911388, 1.7786765197552974,
                            1.7932827571744239, 1.8144097528721321, 1.8374030647656392,
                            1.8549480116192634, 1.8843828007032664, 1.9022522283685337,
                            1.9252154654977307, 1.943266670499853, 1.9530994128293964,
                            1.9523395631811735, 1.9388867623580621, 1.8920765906408807,
                            1.7862576901644667]
                            )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                    [3.849206405998183, 3.9766491772396315, 4.034722228432456,
                    4.097461249180057, 4.12339763790851, 4.145641357994914,
                    4.158394587035046, 4.162759818312447, 4.162640308461358,
                    4.158652791540587, 4.150244798448067, 4.13843460359423,
                    4.128154119612607, 4.108442742835804, 4.095053497479504,
                    4.075550146384204, 4.054802877889783, 4.038745427314903,
                    4.012462183002898, 3.995904929903222, 3.9734860500800373,
                    3.951121361487323, 3.934618210051378, 3.90738706156939,
                    3.8897138188443643, 3.863953151383705, 3.835490695616737,
                    3.8130642834181403, 3.777671522763649, 3.758127431087184,
                    3.7365555371395347, 3.721120808219133, 3.7137189420997623,
                    3.708438968317632, 3.709237851993837, 3.7152117935334936,
                    3.7272017741964327, 3.740011010946214, 3.7676861912386155,
                    3.788398577738892, 3.818689458877156, 3.8484607182632082,
                    3.868878440852897, 3.8987387053259166, 3.91606188756669,
                    3.9386883290220798, 3.9610389199846865, 3.977551992365722,
                    4.004293748816329, 4.020733047913768, 4.042788444640354,
                    4.064157418032696, 4.079290010120359, 4.10188212549747,
                    4.115007853537036, 4.130714757666521, 4.144038649347087,
                    4.151882461843888, 4.161008730272541, 4.1634440412857,
                    4.162174065725458, 4.153752488168884, 4.142684635483829,
                    4.111500095022257, 4.081493913053946, 4.021130531360668,
                    3.935348471939645]
                   )
        end
    end

    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end

using .coll_krook_tests

coll_krook_tests.runtests()
