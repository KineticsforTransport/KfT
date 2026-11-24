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
                 [-0.1902124253234512, -0.03796995693017286, -0.00570107099343619,
                  0.02808895522376095, 0.03440212509803101, 0.028116250580908733,
                  0.023360627125581172, 0.018339365646373, 0.009103782847940263,
                  0.003817111272748658, -0.003736292176186553, -0.010651607956570988,
                  -0.014680323194221398, -0.02039402897174744, -0.022453564823288093,
                  -0.02406143383450459, -0.023913867050473357, -0.022900332528569198,
                  -0.0193914245600555, -0.016745724695217578, -0.01211080968840595,
                  -0.006912064213632023, -0.0032391855836581534, 0.003543367322899518,
                  0.007152554168671884, 0.01225395530118123, 0.016942993196192965,
                  0.019813673110162806, 0.02438402969882035, 0.02646828352773166,
                  0.029004901665822744, 0.030852451045864607, 0.03169521885017839,
                  0.03238079434569957, 0.032287864497472, 0.031540380014113044,
                  0.03006792868191331, 0.02867672937357503, 0.02532324577946861,
                  0.023072016387804966, 0.019214047109921695, 0.014787291663382758,
                  0.011533299957399778, 0.00511919444425127, 0.0014596699991376455,
                  -0.004040821955828618, -0.009417792554408269, -0.012830386842952748,
                  -0.018300676038965204, -0.020680129785596805, -0.023177052267863346,
                  -0.0241960755780187, -0.023933930944787664, -0.021381497058251347,
                  -0.018876737180178814, -0.013865222613716983, -0.007489507300092801,
                  -0.0026856456146849934, 0.006802991683316314, 0.012062246872276049,
                  0.019415197556208172, 0.026277089949185307, 0.02888050644422542,
                  0.032784010087357776, 0.020020289286785203, -0.011887226338711273,
                  -0.07096228071585224]
                )
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                 [1.582947079333887, 1.849370197018268, 1.9057204425486938,
                  1.9481716501216166, 1.9567693877656127, 1.9516743305271995,
                  1.9368234637669894, 1.9225061291827377, 1.8933420179886349,
                  1.8760777773324762, 1.8512082826920981, 1.826945511198538,
                  1.8116286991400523, 1.7852510215714692, 1.772255664505168,
                  1.7543259321108455, 1.7384531523536975, 1.7286394777282916,
                  1.712894472952926, 1.7052981976717878, 1.695606522572935,
                  1.6875118892637506, 1.68304057987595, 1.6761379516286095,
                  1.6728608185584237, 1.668097375473337, 1.6639049862860578,
                  1.6612843597986247, 1.657277932059456, 1.6553720404989953,
                  1.6531840711036017, 1.6515033982218754, 1.6508275118923477,
                  1.6501768339967626, 1.6503142156842014, 1.6509043849243032,
                  1.6522479283342006, 1.6534157853413334, 1.656445189802225,
                  1.6583832813464843, 1.6618770089046575, 1.6657884829280232,
                  1.6688188268103978, 1.6746725329160936, 1.6781809449444833,
                  1.6839038880462842, 1.6911528653558001, 1.6968626427111624,
                  1.7095302178086813, 1.7174444644403775, 1.7307520908673313,
                  1.7456760245482812, 1.7568989509402533, 1.779483886614038,
                  1.7931774203028952, 1.8147940140423917, 1.83832004389535,
                  1.8546261648926272, 1.8859341141108292, 1.9024747099803818,
                  1.9260126267022741, 1.9445349220192296, 1.9539928373535975,
                  1.9531983749799064, 1.9415982940258076, 1.8947106527294562,
                  1.7869062278574823]
                )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                 [3.850876142937179, 3.9811782403748843, 4.035710794245257,
                  4.10175314697682, 4.1256532649027156, 4.148298640533101,
                  4.1608988197124175, 4.164671463330661, 4.164307985477543,
                  4.160418447052964, 4.151735318231826, 4.139502878621313,
                  4.129909189222164, 4.108961738541581, 4.096511620480113,
                  4.076519017137139, 4.055292223810297, 4.040331314987474,
                  4.012469141883595, 3.9970322498461046, 3.9741432076886802,
                  3.951297427274351, 3.935936629792068, 3.9070753593051877,
                  3.890641598856295, 3.8643815558558376, 3.8353587935019107,
                  3.8145297188623313, 3.7771281467968567, 3.7589687155753775,
                  3.7369289759275364, 3.721269607222204, 3.7143697538122247,
                  3.708752888234208, 3.7095449313900013, 3.715598233300949,
                  3.7278743418049363, 3.7397206432525416, 3.769002916112195,
                  3.788269009961944, 3.819103667364665, 3.8493944586640048,
                  3.868399635601561, 3.900033596155913, 3.9161485047519564,
                  3.9392527730980618, 3.962088432095664, 3.977455033247807,
                  4.005807689666145, 4.0211060140295345, 4.0436275647364,
                  4.06542215700367, 4.079588139373767, 4.103610872361545,
                  4.1158778875157465, 4.131997326492202, 4.1457322890464186,
                  4.153089991259052, 4.162929587340712, 4.165197662356767,
                  4.164238298898817, 4.155967294740156, 4.14595996389889,
                  4.113066457337312, 4.085653514646325, 4.024462721235572,
                  3.936174362593473]
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
