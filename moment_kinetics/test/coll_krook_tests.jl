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
                 [-0.19037417829879977, -0.0381143341215827, -0.005835871722714732,
                  0.027970609108040296, 0.03431511883265422, 0.028026675616975034,
                  0.023276006453438004, 0.01825282075926404, 0.009018587700762762,
                  0.003731504807565948, -0.0038219339581452684, -0.01073783460522744,
                  -0.014767121294815133, -0.020482080347597742, -0.022542360061183202,
                  -0.024151231096428444, -0.024004488679279075, -0.022991387311364575,
                  -0.01948304103014207, -0.01683748861283654, -0.012202567629141547,
                  -0.007003609980182551, -0.003330524758255823, 0.0034525057159478277,
                  0.0070619976734974475, 0.012163869824008893, 0.016853336712715963,
                  0.019724236191019838, 0.024294920804910938, 0.02637933301090776,
                  0.028916156368812902, 0.030763813730486922, 0.031606590161855755,
                  0.03229215728360314, 0.0321992287390142, 0.031451751747002425,
                  0.029979257390235298, 0.028587957757616717, 0.02523420734059935,
                  0.022982810915670653, 0.019124567008803925, 0.014697445697149918,
                  0.011443145833069492, 0.005028462727172737, 0.0013686468468864631,
                  -0.004132209100736628, -0.0095094570561638, -0.012922161270532712,
                  -0.018392371735678117, -0.02077160890795592, -0.02326802074490072,
                  -0.024286332588730936, -0.02402359337068816, -0.02146986198416997,
                  -0.018964412758298875, -0.013951866752448833, -0.007575511410812568,
                  -0.002771186029118007, 0.00671742405784485, 0.011977220873096975,
                  0.019329107422766043, 0.02618907083788745, 0.028793326548862398,
                  0.03267687325815859, 0.019893897045017442, -0.01202439978839079,
                  -0.07111645528443215]
                )
        end
        @testset "coll_krook_test_n10" begin
            dkions_n10["output"]["base_directory"] = test_output_directory
            coll_krook_n10["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n10, coll_krook_n10,
                 [1.5828659693714382, 1.8493115807392446, 1.905666677077225,
                  1.9481220217069275, 1.9567205350849353, 1.951623979909503,
                  1.9367703229017463, 1.9224507142328737, 1.8932829808397984,
                  1.876017025440367, 1.851145675501519, 1.8268815418344406,
                  1.8115641686142667, 1.7851857723533557, 1.772190256694428,
                  1.754260359892857, 1.7383876396533844, 1.7285739354131184,
                  1.7128289651226212, 1.70523270282826, 1.6955409447023908,
                  1.6874461135999717, 1.6829747403222868, 1.6760718855842622,
                  1.6727943407275676, 1.6680295739018725, 1.663836653503799,
                  1.6612164702993875, 1.65721118846505, 1.6553053984708557,
                  1.6531173109226014, 1.6514365999320082, 1.6507607713134178,
                  1.6501099787251992, 1.6502473985895412, 1.6508376084887075,
                  1.6521811624414835, 1.6533489990917047, 1.6563785424344077,
                  1.6583162676330296, 1.6618090171132112, 1.6657201590616992,
                  1.6687512546313803, 1.6746063266377331, 1.678115035254357,
                  1.683838011384358, 1.6910872696665196, 1.6967970313425333,
                  1.709464716797191, 1.717378952801199, 1.7306865802861389,
                  1.7456104554913072, 1.7568334251755005, 1.7794185382372505,
                  1.7931123725254843, 1.81472955476512, 1.8382566746265216,
                  1.854563745352071, 1.8858743172407129, 1.9024166644968075,
                  1.9259577566841988, 1.944483146251994, 1.9539429636856647,
                  1.953149248709611, 1.9415480280085902, 1.894655864476657,
                  1.7868440767001754]
                )
        end
        @testset "coll_krook_test_n100" begin
            dkions_n100["output"]["base_directory"] = test_output_directory
            coll_krook_n100["output"]["base_directory"] = test_output_directory
            run_test_with_restart(dkions_n100, coll_krook_n100,
                 [3.8508759993349537, 3.981177762457864, 4.035710066502857,
                  4.101750107108576, 4.125649556386772, 4.14829415225645,
                  4.160894087229428, 4.164666755675208, 4.164303468820065,
                  4.160414030300517, 4.151731148378273, 4.139498886403313,
                  4.129905305739464, 4.108957793685778, 4.096507618793927,
                  4.07651473295869, 4.055287556323361, 4.0403262646176215,
                  4.012463400384953, 3.9970259851897048, 3.9741359110163827,
                  3.9512884835799436, 3.9359260471115607, 3.907059616750283,
                  3.8906214431239095, 3.864352404235802, 3.835320127216834,
                  3.814487034935497, 3.7770873823822484, 3.75893279537091,
                  3.7369014039197115, 3.721248716291318, 3.7143518782606244,
                  3.7087373214064474, 3.7095290820870797, 3.7155797928966283,
                  3.727850625267349, 3.7396918795978573, 3.768964065462492,
                  3.788226416689664, 3.819061599844085, 3.8493600239731123,
                  3.8683719526440465, 3.9000160896138056, 3.916134737403999,
                  3.9392425798177144, 3.962080372965939, 3.9774479025731235,
                  4.005801754040035, 4.021100492153663, 4.043622620466577,
                  4.065417668523332, 4.079583925986081, 4.103606886638406,
                  4.115874011244611, 4.131993395572134, 4.145728233321079,
                  4.153085760471984, 4.162925117866343, 4.16519306454267,
                  4.164233590398517, 4.15596262208178, 4.145955564890522,
                  4.113063109422062, 4.085651165087062, 4.024462066917823,
                  3.9361744681219144]
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
