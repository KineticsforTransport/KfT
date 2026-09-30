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
    test_output_directory = get_MPI_tempdir()

    @testset "coll_krook tests" verbose=use_verbose begin
        println("coll_krook tests")
        @testset "coll_krook_test_n1.0" begin
            dkions_n1["output"]["base_directory"] = test_output_directory
            coll_krook_n1["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n1, coll_krook_n1,
                        [-0.2046469802407029, -0.05530425831794206, -0.020117501842927252,
                         0.014079635383329745, 0.02353343347504934, 0.020308533125974505,
                         0.016885304722402187, 0.012863146834982365,
                         0.0058363528915800055, 0.0012113982148562602,
                         -0.0047902013442442486, -0.01018420104167577,
                         -0.013477464752367283, -0.01712738416714984,
                         -0.01821990182096412, -0.018192400695025934,
                         -0.01637879450955101, -0.013965713875611636,
                         -0.00841343122337383, -0.004120420255472046,
                         0.002462142185913498, 0.00963082666758919, 0.015125563058960113,
                         0.023991652680890246, 0.02935301300548648, 0.03626143260055453,
                         0.04263151759829557, 0.04691632563850952, 0.05293439881117573,
                         0.05607369316429061, 0.05954898736580089, 0.06209903064032829,
                         0.06336936274972892, 0.06428491070406663, 0.0641488258219465,
                         0.06310430922888458, 0.06108330625703043, 0.05898143892770187,
                         0.05454683432281046, 0.051171196788727336, 0.045893656005686824,
                         0.03987867211193884, 0.03503554774168609, 0.026665268719433217,
                         0.021229990844868714, 0.013758382993843613,
                         0.0063958983155080265, 0.0012046760625744563,
                         -0.006360362428714957, -0.010352778355924467,
                         -0.014640998030679806, -0.017397964836854737,
                         -0.018327542489476458, -0.017783789973848985,
                         -0.01621914028020716, -0.012738960630284586,
                         -0.007890345154809774, -0.0037536060220340166,
                         0.0035434806444787307, 0.008130570885893362,
                         0.013886033213720042, 0.018929071112890907, 0.021006487713828546,
                         0.02087616550898478, 0.005205697435021778, -0.027840842777426007,
                         -0.08764314813564641]
                       )
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                        [0.8516879662763227, 1.1162779545555745, 1.1580549327312921,
                         1.1803401815512058, 1.1853213208894406, 1.180195858519735,
                         1.1716568978439759, 1.1625484237660053, 1.1467616275950705,
                         1.1360995605194806, 1.1229702078033323, 1.1105545035580686,
                         1.1028330433426068, 1.092704902035473, 1.0867472261335078,
                         1.0810536787492362, 1.0799264308202625, 1.082414956186776,
                         1.0930278230458168, 1.1034423275247034, 1.1215465752817018,
                         1.1427588244873543, 1.1596312091827112, 1.187274533460897,
                         1.204049967096477, 1.225475777797852, 1.2449770372942135,
                         1.2578665889872682, 1.275709059827298, 1.2848522788205767,
                         1.2948588797291085, 1.3020683217132172, 1.305633229401681,
                         1.3081815396782401, 1.3078084731385156, 1.3048798878750687,
                         1.2992114068169696, 1.2932242272705061, 1.2804231932157908,
                         1.270512685160619, 1.2548187824554966, 1.2365920081086983,
                         1.2217088786353, 1.1956430040677366, 1.1786458355175058,
                         1.1553773788811845, 1.1330500312359593, 1.1179125616855379,
                         1.0978441382881825, 1.0888046330326018, 1.0815846368078457,
                         1.0798060934764444, 1.08189749970281, 1.0896820478894316,
                         1.09583421956809, 1.1044448417467345, 1.1158884570313485,
                         1.1250102060124216, 1.1414990324315812, 1.1516590328134109,
                         1.1654097420929785, 1.1758861807161682, 1.1827851824323634,
                         1.1837397579390092, 1.1775353587409154, 1.1478930962735885,
                         1.074454073538716]
                       )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                    [2.9972006605553587, 3.2783464962839894, 3.329048881435268,
                     3.351713863560354, 3.347795586460847, 3.3271884040340667,
                     3.3082089082544606, 3.2928597284031498, 3.2693793430010305,
                     3.255945836710597, 3.244111480143473, 3.2396471206377875,
                     3.242695330236158, 3.2591019951966644, 3.2763906699476553,
                     3.306401692740126, 3.3420967531560404, 3.370409074440855,
                     3.416912044415177, 3.4449658280806363, 3.480785501344436,
                     3.513206339279647, 3.5347215702351455, 3.5648072814646654,
                     3.580723870750244, 3.5990909560743685, 3.6141486420894022,
                     3.6233547828855515, 3.6351834329605333, 3.640867906277447,
                     3.6467958622942636, 3.6508977994430825, 3.6528661785102226,
                     3.654257065588487, 3.654052037273477, 3.6524581374057976,
                     3.6492893188201507, 3.6458523720839615, 3.638144030685176,
                     3.6318450520176735, 3.62122251009127, 3.607851172542727,
                     3.5959973332839983, 3.572937869622078, 3.5559782657950247,
                     3.5295898159018755, 3.499272587407631, 3.4744503961141735,
                     3.4309310870323872, 3.40238082514316, 3.3633514630079695,
                     3.3257336193170426, 3.3004952529257516, 3.2670364648631054,
                     3.252472205258913, 3.2412889095897244, 3.2405755543055004,
                     3.24537240611741, 3.2623788612189566, 3.2761501240751194,
                     3.297252534880949, 3.3164802844422567, 3.33284227558094,
                     3.351992119298723, 3.351857300783173, 3.3164960033268187,
                     3.2210757251280544]
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
