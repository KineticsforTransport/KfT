module KineticElectronsTests

# Regression test with kinetic electrons, using wall boundary conditions, with recycling
# fraction less than 1 and a plasma source. Runs a Boltzmann electron simulation, restarts
# as a kinetic electron simulation, and checks the final Ez profile.

include("setup.jl")

using moment_kinetics.communication
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info,
                                 postproc_load_variable
using moment_kinetics.looping

using moment_kinetics.Glob

# Input for Boltzmann electron part of run
boltzmann_input = OptionsDict(
    "output" => OptionsDict("run_name" => "kinetic_electron_test_boltzmann_initialisation"),
    "evolve_moments" => OptionsDict("density" => true,
                                    "parallel_flow" => true,
                                    "pressure" => true),
    "r" => OptionsDict("ngrid" => 1,
                       "nelement" => 1),
    "z" => OptionsDict("ngrid" => 5,
                       "discretization" => "gausslegendre_pseudospectral",
                       "nelement" => 8,
                       "bc" => "wall"),
    "vpa" => OptionsDict("ngrid" => 6,
                         "discretization" => "gausslegendre_pseudospectral",
                         "nelement" => 17,
                         "L" => 41.569219381653056,
                         "bc" => "zero",
                         "element_spacing_option" => "coarse_tails8.660254037844386"),
    "composition" => OptionsDict("T_e" => 0.2,
                                 "n_ion_species" => 1,
                                 "n_neutral_species" => 0),
    "ion_species_1" => OptionsDict("initial_temperature" => 0.06666666666666667,
                                   "initial_density" => 1.0),
    "z_IC_ion_species_1" => OptionsDict("initialization_option" => "gaussian",
                                        "density_amplitude" => 1.0,
                                        "temperature_amplitude" => 0.0,
                                        "density_phase" => 0.0,
                                        "upar_amplitude" => 1.4142135623730951,
                                        "temperature_phase" => 0.0,
                                        "upar_phase" => 0.0),
    "vpa_IC_ion_species_1" => OptionsDict("initialization_option" => "gaussian",
                                          "density_amplitude" => 1.0,
                                          "temperature_amplitude" => 0.0,
                                          "density_phase" => 0.0,
                                          "upar_amplitude" => 0.0,
                                          "temperature_phase" => 0.0,
                                          "upar_phase" => 0.0),
    "krook_collisions" => OptionsDict("use_krook" => true),
    "reactions" => OptionsDict("electron_ionization_frequency" => 0.0,
                               "ionization_frequency" => 0.7071067811865476,
                               "charge_exchange_frequency" => 1.0606601717798214),
    "ion_source_1" => OptionsDict("active" => true,
                                  "z_profile" => "gaussian",
                                  "z_width" => 0.25,
                                  "source_strength" => 2.8284271247461903,
                                  "source_T" => 2.0),
    "ion_source_2" => OptionsDict("active" => true,
                                  "z_profile" => "wall_exp_decay",
                                  "z_width" => 0.25,
                                  "source_strength" => 0.7071067811865476,
                                  "source_T" => 0.2),
    "timestepping" => OptionsDict("type" => "SSPRK4",
                                  "nstep" => 20000,
                                  "dt" => 7.071067811865475e-5,
                                  "nwrite" => 2500,
                                  "nwrite_dfns" => 2500,
                                  "steady_state_residual" => true),
    "ion_numerical_dissipation" => OptionsDict("force_minimum_pdf_value" => 0.0),
    "electron_numerical_dissipation" => OptionsDict("force_minimum_pdf_value" => 0.0),
   )

# Test use distributed-memory when possible
if global_size[] % 2 == 0
    # Divide by 2 so that we use shared memory when running in parallel, and so test the
    # ADI preconditioner.
    procs_to_divide_by = global_size[] ÷ 2
else
    procs_to_divide_by = global_size[]
end
boltzmann_input["z"]["nelement_local"] = boltzmann_input["z"]["nelement"] ÷ gcd(boltzmann_input["z"]["nelement"], procs_to_divide_by)

kinetic_input = deepcopy(boltzmann_input)
kinetic_input["output"]["run_name"] = "kinetic_electron_test"
kinetic_input["composition"]["electron_physics"] = "kinetic_electrons"
kinetic_input["timestepping"] = OptionsDict("type" => "PareschiRusso2(2,2,2)",
                                            "kinetic_electron_solver" => "implicit_p_implicit_pseudotimestep",
                                            "kinetic_electron_preconditioner" => "lu",
                                            "kinetic_ion_solver" => "full_explicit_ion_advance",
                                            "nstep" => 100,
                                            "dt" => 7.0710678118654756e-6,
                                            "nwrite" => 100,
                                            "nwrite_dfns" => 100,
                                           )

kinetic_input["electron_timestepping"] = OptionsDict("nstep" => 5000000,
                                                     "dt" => 1.4142135623730951e-5,
                                                     "maximum_dt" => Inf,
                                                     "nwrite" => 10000,
                                                     "nwrite_dfns" => 100000,
                                                     "decrease_dt_iteration_threshold" => 5000,
                                                     "increase_dt_iteration_threshold" => 0,
                                                     "cap_factor_ion_dt" => 10.0,
                                                     "initialization_residual_value" => 1.0e0,
                                                     "converged_residual_value" => 0.1,
                                                    )

kinetic_input["nonlinear_solver"] = OptionsDict("nonlinear_max_iterations" => 1000,
                                                "rtol" => 1.0e-8,
                                                "atol" => 1.0e-14,
                                                "linear_restart" => 5,
                                                "preconditioner_update_interval" => 100,
                                               )

kinetic_input_adaptive_timestep = deepcopy(kinetic_input)
kinetic_input_adaptive_timestep["output"]["run_name"] = "kinetic_electron_adaptive_timestep_test"
kinetic_input_adaptive_timestep["timestepping"]["type"] = "KennedyCarpenterARK324"
kinetic_input_adaptive_timestep["timestepping"]["maximum_dt"] = 7.0710678118654756e-6

kinetic_input_time_evolving = deepcopy(kinetic_input_adaptive_timestep)
kinetic_input_time_evolving["output"]["run_name"] = "kinetic_electron_implicit_time_evolving_test"
kinetic_input_time_evolving["timestepping"]["kinetic_electron_solver"] = "implicit_time_evolving"

# The following settings don't run for long enough to give a very good test, but at least
# make sure that the explicit electron solver can take a single timestep.
kinetic_input_implicit_ppar_explicit_pseudotimestep = deepcopy(kinetic_input)
kinetic_input_implicit_ppar_explicit_pseudotimestep["output"]["run_name"] = "kinetic_electron_implicit_ppar_explicit_pseudotimestep_test"
kinetic_input_implicit_ppar_explicit_pseudotimestep["timestepping"]["nstep"] = 1
kinetic_input_implicit_ppar_explicit_pseudotimestep["timestepping"]["kinetic_electron_solver"] = "implicit_p_explicit_pseudotimestep"
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["nstep"] = 5000000
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["dt"] = 7.071067811865474e-8
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["maximum_dt"] = 7.0710678118654756e-6
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["nwrite"] = 10000
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["nwrite_dfns"] = 100000
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["type"] = "Fekete4(3)"
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["rtol"] = 0.001
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["atol"] = 1.0e-14
kinetic_input_implicit_ppar_explicit_pseudotimestep["electron_timestepping"]["minimum_dt"] = 7.071067811865476e-10

kinetic_input_explicit_time_evolving = deepcopy(kinetic_input_implicit_ppar_explicit_pseudotimestep)
kinetic_input_explicit_time_evolving["output"]["run_name"] = "kinetic_electron_explicit_time_evolving_test"
kinetic_input_explicit_time_evolving["timestepping"]["kinetic_electron_solver"] = "explicit_time_evolving"
kinetic_input_explicit_time_evolving["timestepping"]["type"] = "Fekete4(3)"
kinetic_input_explicit_time_evolving["timestepping"]["maximum_dt"] = 7.0710678118654756e-6

"""
Run a test for a single set of parameters
"""
function run_test()
    test_output_directory = get_MPI_tempdir()

    this_boltzmann_input = deepcopy(boltzmann_input)
    this_boltzmann_input["output"]["base_directory"] = test_output_directory

    # Suppress console output while running.
    quietoutput() do
        run_moment_kinetics(this_boltzmann_input)
    end

    preconditioners = ["schur_complement", "lu"]

    # Test implicit electron solve
    test_inputs = [(deepcopy(kinetic_input), "fixed timestep", 1.0e-6),
                   (deepcopy(kinetic_input_time_evolving), "time evolving", 4.0e-4),]
    @long push!(test_inputs, (deepcopy(kinetic_input_adaptive_timestep), "adaptive timestep", 2.0e-4))
    @testset "$label $precon" for (this_kinetic_input, label, tol) ∈ test_inputs,
                                  precon ∈ preconditioners

        this_kinetic_input["output"]["base_directory"] = test_output_directory
        this_kinetic_input["output"]["run_name"] *= "_$precon"
        this_kinetic_input["timestepping"]["kinetic_electron_preconditioner"] = precon

        println("    - testing kinetic electrons, $label, $precon")
        if precon == "lu"
            # Provide some progress info
            # Don't use distributed memory parallelism with LU preconditioner
            pop!(this_kinetic_input["z"], "nelement_local")
        end

        # Suppress console output while running.
        quietoutput() do
            restart_from_directory = joinpath(this_boltzmann_input["output"]["base_directory"], this_boltzmann_input["output"]["run_name"])
            restart_from_file_pattern = this_boltzmann_input["output"]["run_name"] * ".dfns*.h5"
            restart_from_file = glob(restart_from_file_pattern, restart_from_directory)[1]

            # run kinetic electron simulation
            run_moment_kinetics(this_kinetic_input; restart=restart_from_file)
        end

        if global_rank[] == 0
            # Load and analyse output
            #########################

            path = joinpath(realpath(this_kinetic_input["output"]["base_directory"]), this_kinetic_input["output"]["run_name"])

            # open the output file(s)
            run_info = get_run_info_no_setup(path, dfns=true)

            # load fields data
            Ez = postproc_load_variable(run_info, "Ez")[:,1,:]
            vthe = postproc_load_variable(run_info, "electron_thermal_speed")[:,1,:]
            electron_advance_linear_iterations = postproc_load_variable(run_info, "electron_advance_linear_iterations")[end]

            close_run_info(run_info)

            # Regression test
            # Benchmark data generated in serial on Linux with the LU preconditioner
            expected_Ez = [-1.1780452305101836 -1.1829749301103147;
                           -1.0263898136838296 -1.0347132584217122;
                           -0.6619800344685941 -0.671967087143512;
                           -0.42769948572804395 -0.43779408866897335;
                           -0.43328104167393633 -0.4436714199811528;
                           -0.34258386854929335 -0.3527296800430666;
                           -0.2846611112216824 -0.29406073863004306;
                           -0.2632585878015962 -0.27143843617728214;
                           -0.24375953332618883 -0.2515128828023723;
                           -0.2439294490302325 -0.2507180093002268;
                           -0.22556194856501618 -0.23112630401136958;
                           -0.19601625359487285 -0.2003596763350178;
                           -0.17351491614712533 -0.17736863139793682;
                           -0.15163598429416322 -0.1547354496144346;
                           -0.09745645498055537 -0.09937914326643132;
                           -0.03451420674828861 -0.03518385621268484;
                           2.579009934951453e-14 6.559383356841287e-15;
                           0.03451420674828553 0.03518385621269143;
                           0.09745645498057001 0.09937914326645703;
                           0.1516359842941711 0.15473544961442073;
                           0.17351491614708064 0.1773686313978434;
                           0.19601625359485825 0.2003596763350174;
                           0.22556194856501882 0.2311263040113987;
                           0.24392944903023425 0.2507180093001983;
                           0.2437595333261215 0.2515128828023129;
                           0.26325858780158995 0.27143843617728564;
                           0.2846611112217038 0.2940607386300604;
                           0.34258386854927847 0.3527296800430507;
                           0.43328104167393244 0.44367141998115317;
                           0.4276994857280369 0.4377940886689848;
                           0.6619800344685978 0.6719670871435005;
                           1.0263898136838518 1.0347132584217178;
                           1.1780452305101727 1.1829749301103543]
            expected_vthe = [18.679721002547694 18.537399093466572;
                             19.527570628031825 19.409938211549623;
                             20.710503446243816 20.631423377737665;
                             21.426455865241746 21.376629155996586;
                             21.682612092310634 21.64706315108888;
                             21.90170057132989 21.879324177030913;
                             22.20559380695977 22.206115827571413;
                             22.431812913381506 22.4515226114792;
                             22.523871086117758 22.552361202895167;
                             22.613380694199073 22.649612927234486;
                             22.7510357713038 22.79909959927407;
                             22.86647049445335 22.9235636934246;
                             22.912885888201252 22.973497625425527;
                             22.958097299199483 23.021985437610567;
                             23.016580611165626 23.084383028531217;
                             23.048846391201412 23.11885428254528;
                             23.050567053530717 23.12056676677472;
                             23.04884639120141 23.118854282545275;
                             23.016580611165626 23.084383028531217;
                             22.95809729919949 23.021985437610557;
                             22.912885888201263 22.973497625425516;
                             22.86647049445335 22.923563693424576;
                             22.75103577130381 22.799099599274037;
                             22.613380694199073 22.64961292723447;
                             22.52387108611776 22.552361202895217;
                             22.431812913381506 22.45152261147921;
                             22.205593806959772 22.20611582757142;
                             21.90170057132989 21.879324177030913;
                             21.68261209231064 21.64706315108884;
                             21.426455865241742 21.376629155996557;
                             20.710503446243813 20.631423377737622;
                             19.527570628031818 19.40993821154964;
                             18.679721002547698 18.537399093466583]

            if expected_Ez == nothing
                # Error: no expected input provided
                println("data tested would be: Ez=", Ez)
                @test false
            else
                @test elementwise_isapprox(Ez, expected_Ez, rtol=0.0, atol=2.0*tol)
            end
            if expected_vthe == nothing
                # Error: no expected input provided
                println("data tested would be: vthe=", vthe)
                @test false
            else
                @test elementwise_isapprox(vthe, expected_vthe, rtol=tol, atol=0.0)
            end

            if this_kinetic_input["timestepping"]["kinetic_electron_preconditioner"] != "adi"
                # Iteration counts are fairly inconsistent, but it's good to check that they at
                # least don't unexpectedly increase by an order of magnitude.
                # Expected iteration count is from a serial run on Linux.
                # ADI preconditioner requires more iterations, so don't check iteration
                # counts when using ADI.
                expected_electron_advance_linear_iterations = 13763
                @test electron_advance_linear_iterations < 2 * expected_electron_advance_linear_iterations
                if !(electron_advance_linear_iterations < 2 * expected_electron_advance_linear_iterations)
                    println("electron_advance_linear_iterations=$electron_advance_linear_iterations was greater than twice the expected $expected_electron_advance_linear_iterations.")
                end
            end
        end
    end

    # Test explicit electron solves - separate loop for different expected results due to
    # a shorter simulation time, and no test for iteration counts.
    @long @testset "$label" for (this_kinetic_input, label, tol) ∈ ((deepcopy(kinetic_input_implicit_ppar_explicit_pseudotimestep), "explicit solve", 1.0e-4),
                                                                    (deepcopy(kinetic_input_explicit_time_evolving), "explicit time evolving", 4.0e-4),)
        this_kinetic_input["output"]["base_directory"] = test_output_directory

        # Provide some progress info
        println("    - testing kinetic electrons, $label")

        # Suppress console output while running.
        quietoutput() do
            restart_from_directory = joinpath(this_boltzmann_input["output"]["base_directory"], this_boltzmann_input["output"]["run_name"])
            restart_from_file_pattern = this_boltzmann_input["output"]["run_name"] * ".dfns*.h5"
            restart_from_file = glob(restart_from_file_pattern, restart_from_directory)[1]

            # run kinetic electron simulation
            run_moment_kinetics(this_kinetic_input; restart=restart_from_file)
        end

        if global_rank[] == 0
            # Load and analyse output
            #########################

            path = joinpath(realpath(this_kinetic_input["output"]["base_directory"]), this_kinetic_input["output"]["run_name"])

            # open the output file(s)
            run_info = get_run_info_no_setup(path, dfns=true)

            # load fields data
            Ez = postproc_load_variable(run_info, "Ez")[:,1,:]
            vthe = postproc_load_variable(run_info, "electron_thermal_speed")[:,1,:]

            close_run_info(run_info)

            # Regression test
            # Benchmark data generated in serial on Linux
            expected_Ez = [-1.1262819815582732 -1.1268876490806572;
                           -1.0185605525334143 -1.0188208753017798;
                           -0.6741556175828979 -0.6741843327070705;
                           -0.43456846812123523 -0.4347115554671706;
                           -0.44851712786767417 -0.44864031342970956;
                           -0.3546417200934376 -0.3546995721979947;
                           -0.2949673971075837 -0.2950668535051838;
                           -0.2720667189401215 -0.272142459735386;
                           -0.25144481472004965 -0.25152779923741053;
                           -0.251157338372241 -0.25122581216340156;
                           -0.23160290173033277 -0.2316589899788601;
                           -0.20084171520879382 -0.20088524042085756;
                           -0.17756128926133466 -0.1776003616392862;
                           -0.1552431430057613 -0.15527381644509414;
                           -0.09973568387974677 -0.09975463885205996;
                           -0.03529768576136534 -0.03530427967923867;
                           3.607632394556368e-14 3.205128203201748e-14;
                           0.03529768576135704 0.03530427967923024;
                           0.09973568387976896 0.09975463885208218;
                           0.15524314300577274 0.1552738164451057;
                           0.1775612892612843 0.17760036163924134;
                           0.20084171520877006 0.2008852404208342;
                           0.23160290173033965 0.23165898997886872;
                           0.251157338372242 0.2512258121634013;
                           0.25144481471999097 0.25152779923734137;
                           0.272066718940108 0.27214245973537365;
                           0.2949673971076046 0.295066853505208;
                           0.354641720093429 0.3546995721979838;
                           0.44851712786765885 0.4486403134296921;
                           0.43456846812123157 0.4347115554671649;
                           0.6741556175828951 0.6741843327070681;
                           1.0185605525334376 1.0188208753018038;
                           1.1262819815581997 1.1268876490805784]
            expected_vthe = [18.57849595411752 18.57664898782716;
                             19.4023915724272 19.401290309352586;
                             20.618872475200234 20.618172655738483;
                             21.363766968835392 21.3632717274725;
                             21.631857230529953 21.63163758621151;
                             21.870770766810537 21.870573531703418;
                             22.199183505895927 22.19919663340091;
                             22.448202498745466 22.448404790855253;
                             22.54691556670904 22.547204322185078;
                             22.64710133142786 22.647468677479754;
                             22.79644655706123 22.79693257392232;
                             22.923258507986993 22.92383533928831;
                             22.971638057815298 22.972249218669496;
                             23.022540506442624 23.023185205914757;
                             23.084894104913076 23.0855762798366;
                             23.12085244462928 23.121556824992098;
                             23.120880257374072 23.121583234412643;
                             23.120852444629275 23.121556824992098;
                             23.084894104913083 23.08557627983661;
                             23.022540506442606 23.02318520591474;
                             22.971638057815305 22.972249218669504;
                             22.923258507986976 22.923835339288296;
                             22.79644655706124 22.796932573922327;
                             22.64710133142784 22.64746867747973;
                             22.54691556670906 22.547204322185102;
                             22.44820249874545 22.44840479085523;
                             22.199183505895945 22.199196633400927;
                             21.870770766810512 21.87057353170339;
                             21.631857230529977 21.631637586211536;
                             21.363766968835368 21.36327172747248;
                             20.618872475200234 20.61817265573849;
                             19.402391572427184 19.40129030935257;
                             18.578495954117546 18.576648987827188]

            if expected_Ez == nothing
                # Error: no expected input provided
                println("data tested would be: Ez=", Ez)
                @test false
            else
                @test elementwise_isapprox(Ez, expected_Ez, rtol=0.0, atol=2.0*tol)
            end
            if expected_vthe == nothing
                # Error: no expected input provided
                println("data tested would be: vthe=", vthe)
                @test false
            else
                @test elementwise_isapprox(vthe, expected_vthe, rtol=tol, atol=0.0)
            end
        end
    end

    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

function runtests()
    @testset "kinetic electrons" verbose=use_verbose begin
        println("Kinetic electron tests")
        run_test()
    end
    return nothing
end

end # KineticElectronsTests


using .KineticElectronsTests

KineticElectronsTests.runtests()
