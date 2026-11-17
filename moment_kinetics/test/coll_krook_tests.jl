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
     "discretization" => "chebyshev_pseudospectral",
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
                 [-0.19021520855696608, -0.037970195857344555, -0.00570146501664363,
                  0.02808806687994785, 0.034400524392729795, 0.028114787092880275,
                  0.02335944469551856, 0.018338229177454233, 0.009102920309627973,
                  0.003816378556947854, -0.0037369097666849374, -0.010652107061044611,
                  -0.014680737039266277, -0.020394326888769646, -0.022453819757334334,
                  -0.02406161775507584, -0.02391399673937708, -0.022900430532965855,
                  -0.01939148149520511, -0.016745762382816085, -0.012110827218512716,
                  -0.00691206684980497, -0.0032391806869090153, 0.0035433814009812633,
                  0.0071525715392519425, 0.012253975894501723, 0.016943015656697227,
                  0.01981369620539196, 0.024384053372406785, 0.026468307155159212,
                  0.02900492512259144, 0.030852474231785634, 0.03169524194535132,
                  0.03238081732289146, 0.032287887522510744, 0.031540403094673056,
                  0.030067952010440963, 0.02867675282856177, 0.025323269470867084,
                  0.02307203993498566, 0.019214070127558592, 0.014787313364100245,
                  0.011533320223969894, 0.005119210040739532, 0.0014596818486286663,
                  -0.004040818592517518, -0.009417801518894292, -0.012830407059728844,
                  -0.0183007239789662, -0.02068019896263033, -0.023177156497484825,
                  -0.024196228602293472, -0.023934125080326524, -0.0213817733550168,
                  -0.018877071504282195, -0.013865651392215423, -0.00749007594253821,
                  -0.002686275093089654, 0.0068021848699025556, 0.012061310960797555,
                  0.01941403181217924, 0.026275922690565154, 0.028878915301625063,
                  0.032782907443873015, 0.02001961131059853, -0.011887599869571512,
                  -0.07096214715530867]
                )
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                 [1.5829428609672456, 1.8493684592850395, 1.9057187661858401,
                  1.948169849255173, 1.9567676625440567, 1.9516726016918966,
                  1.9368218360929668, 1.922504569858221, 1.893340624764887,
                  1.8760764463453552, 1.8512070480701153, 1.8269443702212813,
                  1.8116276404677032, 1.7852501112477555, 1.7722548449822864,
                  1.7543252548315653, 1.7384526217531115, 1.7286389711552461,
                  1.712893973103113, 1.7052976553704742, 1.695605814759533,
                  1.6875111458508731, 1.6830398805571922, 1.6761378176359785,
                  1.67286100178881, 1.6680978436201856, 1.663905681329,
                  1.6612851944366345, 1.657278946327763, 1.6553731281940882,
                  1.6531852362548278, 1.6515046139808363, 1.650828752601819,
                  1.6501780976589435, 1.6503154717231923, 1.6509056209150623,
                  1.6522491254137386, 1.6534169396120675, 1.6564462366171515,
                  1.6583842460964058, 1.6618778150122573, 1.6657890726643907,
                  1.6688192605248535, 1.6746725417016373, 1.678180629325858,
                  1.6839031634581443, 1.6911521110203047, 1.6968619681692223,
                  1.709529681998837, 1.717444009332654, 1.730751579515947,
                  1.7456754284191494, 1.7568982553786487, 1.7794830132904835,
                  1.7931764664659722, 1.8147929351154237, 1.8383188581833843,
                  1.8546249114015685, 1.8859327500111442, 1.9024732601312786,
                  1.9260110595370687, 1.9445332281658996, 1.9539911354518666,
                  1.9531966122591158, 1.941596583276105, 1.8947089105669876,
                  1.7869046654607665]
                )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                 [3.850824365251345, 3.9811662236750296, 4.035703788674293,
                  4.101748198373712, 4.125649428084572, 4.148295284641993,
                  4.16089603687649, 4.164668942898611, 4.164305968242765,
                  4.160416446007831, 4.151733571875441, 4.13950124116633,
                  4.129907744115463, 4.108960318431994, 4.0965103551647,
                  4.076517860011044, 4.055291202986864, 4.040330300734374,
                  4.012468245248808, 3.9970313224205456, 3.974142208926275,
                  3.9512962553761195, 3.9359353336241214, 3.9070737793171624,
                  3.890639938416305, 3.8643799215914387, 3.83535745228121,
                  3.814528646254379, 3.7771275699966567, 3.758968298602543,
                  3.736928709082885, 3.721269390081526, 3.7143695708281323,
                  3.7087527030055853, 3.709544758257429, 3.715598033913283,
                  3.7278741143933147, 3.739720350406628, 3.769002424578939,
                  3.788268291219241, 3.819102543600086, 3.8493929448027924,
                  3.8683979901816268, 3.900031966957533, 3.916147014897993,
                  3.93925149703474, 3.962087358602485, 3.9774540338145954,
                  4.0058067918252585, 4.021105064695127, 4.043626567825434,
                  4.065421051692561, 4.079587009876378, 4.103609487605383,
                  4.115876518198461, 4.131995792857918, 4.145730649385057,
                  4.1530881383736356, 4.162927605301036, 4.16519541680173,
                  4.164235779272952, 4.155964219942664, 4.145956690809827,
                  4.113062046898235, 4.085648420627556, 4.024454826182284,
                  3.936155355911855]
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
