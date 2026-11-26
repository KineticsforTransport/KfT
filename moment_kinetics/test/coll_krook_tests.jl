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
    test_output_directory = get_MPI_tempdir()

    @testset "coll_krook tests" verbose=use_verbose begin
        println("coll_krook tests")
        @testset "coll_krook_test_n1.0" begin
            dkions_n1["output"]["base_directory"] = test_output_directory
            coll_krook_n1["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n1, coll_krook_n1,
                 [-0.19021109485734375, -0.03796980072389158, -0.005700879804140672,
                  0.02808933227034412, 0.03440275145671171, 0.028116994008074554,
                  0.02336135956454104, 0.018339977593675883, 0.009104400624755208,
                  0.0038177313396465366, -0.0037357481196307035, -0.010651096249092848,
                  -0.014679824266237913, -0.02039358454513802, -0.02245315903102114,
                  -0.024061063089025954, -0.02391353893284694, -0.022900027537179356,
                  -0.019391167120826495, -0.016745489033321345, -0.012110605630159017,
                  -0.006911887722277563, -0.003239025866324484, 0.003543500645726101,
                  0.007152675117608743, 0.01225406035650021, 0.01694308511374283,
                  0.019813757647739006, 0.024384103182356495, 0.02646835224698312,
                  0.029004964782072942, 0.030852510311328532, 0.03169527640650987,
                  0.03238085059670727, 0.032287920907419575, 0.03154043789073454,
                  0.030067989552008734, 0.02867679320835011, 0.025323317091917305,
                  0.023072092971104213, 0.019214133146594377, 0.014787389445278347,
                  0.011533407147125337, 0.005119322230706939, 0.0014598109723210021,
                  -0.00404065876180928, -0.009417603339918678, -0.012830178462621696,
                  -0.01830042799769821, -0.020679858852361798, -0.02317674226645248,
                  -0.0241957272181214, -0.023933555731299124, -0.021381067433680213,
                  -0.018876280555996938, -0.013864716064352326, -0.007489000084184298,
                  -0.002685084602556379, 0.006803610226764141, 0.01206286628129942,
                  0.019415821358442144, 0.026277889582005304, 0.028881230045703112,
                  0.032784487319516394, 0.020020562060828456, -0.011887036996300449,
                  -0.07096225225948793]
                )
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                 [1.5830007561217174, 1.8493777655010057, 1.905724846325108,
                  1.948175766995193, 1.9567717856258804, 1.9516767956817418,
                  1.9368253286192976, 1.9225080386810338, 1.8933434612310032,
                  1.8760793474369983, 1.8512096621341922, 1.8269469521447086,
                  1.811630060778631, 1.7852524819573905, 1.7722570742823802,
                  1.7543274225851289, 1.738454628736336, 1.7286410519231081,
                  1.7128961265633216, 1.7053001172609585, 1.6956087779613314,
                  1.687514724599967, 1.6830436830931963, 1.6761409983068838,
                  1.6728634882946656, 1.66809963296845, 1.6639068991819663,
                  1.6612861094269353, 1.6572794253453076, 1.655373493796307,
                  1.6531853910110446, 1.6515046926955457, 1.6508287597678337,
                  1.650178114170321, 1.6503154523753392, 1.650905668187888,
                  1.652249218831014, 1.65341715058838, 1.6564466429877651,
                  1.6583848738227058, 1.6618787634121586, 1.6657905687749701,
                  1.6688211006905036, 1.6746754517507263, 1.6781840754524913,
                  1.6839069665893704, 1.6911554545916379, 1.6968648483672568,
                  1.709531961077281, 1.717446112252774, 1.7307536029451454,
                  1.7456775400254891, 1.7569003778038519, 1.7794853456669337,
                  1.7931787985026852, 1.8147954486620013, 1.8383214010784241,
                  1.8546276799357946, 1.8859355799178923, 1.9024764185006824,
                  1.9260143737508262, 1.9445372685940332, 1.9539947874832206,
                  1.9532016358471211, 1.941601551564611, 1.8947164128133325,
                  1.7869137315653338]
                )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                 [3.8508963495625452, 3.981191095285182, 4.03572397448814,
                  4.101764283570572, 4.125664036639995, 4.148308691342073,
                  4.160908140087042, 4.1646803649385005, 4.164316291448201,
                  4.160426636326949, 4.15174319649079, 4.139510542436103,
                  4.129916672966457, 4.108969109422947, 4.096518848757869,
                  4.076526076646057, 4.055299087294844, 4.040338102669743,
                  4.0124756926823, 3.997038733629747, 3.974149696830398,
                  3.9513040300273787, 3.935943303882426, 3.9070821010863597,
                  3.8906481303427833, 3.8643872717221956, 3.835363013608726,
                  3.8145328283253264, 3.777129579996112, 3.758969664012403,
                  3.736929553241813, 3.7212700878625933, 3.71437017818448,
                  3.7087533230449887, 3.709545343564215, 3.715598686170298,
                  3.727874839452255, 3.739721273392245, 3.7690040924976267,
                  3.7882708791799447, 3.8191069957349892, 3.8493994572646244,
                  3.868405502540663, 3.9000402883043273, 3.9161552521515923,
                  3.9392594432472956, 3.962094968660344, 3.9774615196719005,
                  4.005814194791248, 4.021112640513887, 4.04363435451506,
                  4.065429126182958, 4.079595200587323, 4.103618206290591,
                  4.115885264721482, 4.132004877060982, 4.145740022474034,
                  4.153097962604323, 4.162937822288554, 4.165206151137801,
                  4.164247242574027, 4.155977024580068, 4.14597005408052,
                  4.113077429544955, 4.085664904755492, 4.024476132248899,
                  3.9361860223203324]
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
