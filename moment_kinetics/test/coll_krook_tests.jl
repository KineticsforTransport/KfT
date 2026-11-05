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
                        [-0.20464635486282068, -0.0553040963811019, -0.020117131077546306,
                         0.014081081952311243, 0.02353884321529961, 0.020310007674498327,
                         0.01689118308116017, 0.012866480381057064, 0.005839650581325177,
                         0.0012140949325176122, -0.004786779876935845,
                         -0.01018129303707813, -0.013474416696112943,
                         -0.01712449390978572, -0.01821720330410044,
                         -0.018189888724736914, -0.016376457346175738,
                         -0.01396356914041006, -0.008411571531252743,
                         -0.004118728285625513, 0.002463626335168406,
                         0.009632115282183619, 0.01512670164594783, 0.023992566837684123,
                         0.029353795635113505, 0.036262080724569296, 0.042632099688389255,
                         0.04691690985536146, 0.052935017444547816, 0.05607431481462667,
                         0.05954957669015937, 0.06209954328982073, 0.06336981544937301,
                         0.06428530993956724, 0.06414923372863199, 0.06310477562113825,
                         0.061083856391715444, 0.05898203797691676, 0.05454745686271487,
                         0.05117180683761481, 0.04589423640908413, 0.03987927373973874,
                         0.03503621591825678, 0.02666611560871056, 0.02123097507577072,
                         0.01375955767460552, 0.006397276834621074, 0.00120619494404891,
                         -0.006358585380339882, -0.010350832201808147,
                         -0.014638808784283766, -0.01739551350651516,
                         -0.018325026480766633, -0.01778093848916719,
                         -0.01621629196482185, -0.012735872345203693,
                         -0.007887531679040403, -0.0037500047892425493,
                         0.003545951833988418, 0.00813471030556842, 0.013889728094449483,
                         0.01893438461552181, 0.02100775792555499, 0.020880120262790763,
                         0.005205752295040232, -0.027840711821963827,
                         -0.08764284946908006]
                       )
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                        [0.8519342339536957, 1.1166115038523083, 1.1584073159788146,
                         1.1806950730706471, 1.185674069140162, 1.1805388460025934,
                         1.1719912389934068, 1.162874812781358, 1.1470767177957484,
                         1.1364066486929179, 1.1232676958228212, 1.1108412183680836,
                         1.1031111551568133, 1.0929683287033267, 1.0870007609537875,
                         1.0812799618895281, 1.0800955816703375, 1.0825198271594558,
                         1.0930478383835216, 1.1034102109754775, 1.1214442444419512,
                         1.1425989443786464, 1.159439289256536, 1.1870428483664532,
                         1.2038002822348586, 1.2252081314269212, 1.2446970715247279,
                         1.257579906733825, 1.27541494303152, 1.2845549271558738,
                         1.2945585865266631, 1.301766034390212, 1.3053301355231561,
                         1.3078778128989414, 1.3075048782276717, 1.3045768798338195,
                         1.2989098912807056, 1.2929243041294454, 1.280127383665726,
                         1.2702205077518305, 1.2545336869036379, 1.236316951992559,
                         1.2214441634153157, 1.1954017719019185, 1.1784253030358174,
                         1.1551925444790425, 1.132913712994914, 1.117822027479836,
                         1.0978394530640847, 1.0888478258825411, 1.0817065079684443,
                         1.0800045112136683, 1.0821313954883605, 1.08994097726845,
                         1.096101887676485, 1.1047247616420897, 1.1161798103654428,
                         1.1253088808720013, 1.141810358193608, 1.1519772993547777,
                         1.1657389913733507, 1.1762243597148867, 1.1831314093117855,
                         1.1840937798111328, 1.1778921254927364, 1.1482402059720298,
                         1.0747704956067758]
                       )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                    [2.9971801507011393, 3.278327705110373, 3.329026381028684,
                     3.3516849052092796, 3.3477636468712597, 3.327148537582419,
                     3.308150341783696, 3.2927791217574303, 3.2692694904943913,
                     3.255828101236384, 3.2439920186713422, 3.2395352973548572,
                     3.2425945478849316, 3.2590244215863344, 3.2763277254455896,
                     3.3063563518425014, 3.342065205642163, 3.370385365391053,
                     3.4168973294100504, 3.44495499507877, 3.4807783634127407,
                     3.5132016771314056, 3.5347181627467843, 3.56480524222563,
                     3.580722384355494, 3.5990899993470347, 3.614148023517342,
                     3.6233543420057597, 3.6351831787868805, 3.6408677311929183,
                     3.6467957562029185, 3.6508977400010396, 3.652866137594857,
                     3.6542570393177614, 3.6540520079378314, 3.6524580947030474,
                     3.6492892413860547, 3.645852257872442, 3.6381438179363976,
                     3.6318447503878653, 3.621222026509097, 3.6078504216793954,
                     3.595996290420223, 3.572936127609955, 3.5559758656061824,
                     3.5295861431232134, 3.499266950591367, 3.4744427051124784,
                     3.4309184270107487, 3.402363713070608, 3.3633259631697445,
                     3.3256964320721063, 3.30044691212718, 3.2669662541483135,
                     3.2523870260137984, 3.241185107566552, 3.2404591777125633,
                     3.2452525899388287, 3.26226420106995, 3.27604712386542,
                     3.2971776779210282, 3.316431841519527, 3.332804498365245,
                     3.351961699288368, 3.3518302232909516, 3.316474341450382,
                     3.2210584927996617]
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
