module BraginskiiElectronsIMEX

# Regression test using wall boundary conditions, with recycling fraction less than 1 and
# a plasma source. Runs for a while and then checks phi profile against saved reference
# output.

include("setup.jl")

using Base.Filesystem: tempname
using MPI

using moment_kinetics.interpolation: interpolate_to_grid_z
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info, get_variable
using moment_kinetics.utils: merge_dict_with_kwargs!

# default inputs for tests
test_input = OptionsDict("composition" => OptionsDict("n_ion_species" => 1,

                                                      "n_neutral_species" => 1,
                                                      "electron_physics" => "braginskii_fluid",
                                                      "T_e" => 0.2),
                         "output" => OptionsDict("run_name" => "braginskii-electrons-imex"),
                         "evolve_moments" => OptionsDict("density" => true,
                                                         "parallel_flow" => true,
                                                         "pressure" => true),
                         "ion_species_1" => OptionsDict("initial_density" => 1.0,
                                                        "initial_temperature" => 0.3333333333333333),
                         "z_IC_ion_species_1" => OptionsDict("initialization_option" => "sinusoid",
                                                             "density_amplitude" => 0.1,
                                                             "density_phase" => 0.0,
                                                             "upar_amplitude" => 1.4142135623730951,
                                                             "upar_phase" => 0.0,
                                                             "temperature_amplitude" => 0.1,
                                                             "temperature_phase" => 1.0),
                         "vpa_IC_ion_species_1" => OptionsDict("initialization_option" => "gaussian",
                                                               "density_amplitude" => 1.0,
                                                               "density_phase" => 0.0,
                                                               "upar_amplitude" => 0.0,
                                                               "upar_phase" => 0.0,
                                                               "temperature_amplitude" => 0.0,
                                                               "temperature_phase" => 0.0),
                         "neutral_species_1" => OptionsDict("initial_density" => 1.0,
                                                            "initial_temperature" => 0.3333333333333333),
                         "z_IC_neutral_species_1" => OptionsDict("initialization_option" => "sinusoid",
                                                                 "density_amplitude" => 0.001,
                                                                 "density_phase" => 0.0,
                                                                 "upar_amplitude" => 0.0,
                                                                 "upar_phase" => 0.0,
                                                                 "temperature_amplitude" => 0.0,
                                                                 "temperature_phase" => 0.0),
                         "vz_IC_neutral_species_1" => OptionsDict("initialization_option" => "gaussian",
                                                                  "density_amplitude" => 1.0,
                                                                  "density_phase" => 0.0,
                                                                  "upar_amplitude" => 0.0,
                                                                  "upar_phase" => 0.0,
                                                                  "temperature_amplitude" => 0.0,
                                                                  "temperature_phase" => 0.0),
                         "electron_fluid_collisions" => OptionsDict("nu_ei" => 1414.213562373095),
                         "reactions" => OptionsDict("charge_exchange_frequency" => 1.0606601717798214,
                                                    "electron_charge_exchange_frequency" => 0.0,
                                                    "ionization_frequency" => 0.7071067811865476,
                                                    "electron_ionization_frequency" => 0.7071067811865476),
                         # Use higher-order "KennedyCarpenterARK437" scheme rather than
                         # "KennedyCarpenterARK324" for the default case in this test,
                         # because the higher order scheme is faster when enforcing the
                         # tight tolerances that we use here. Tight solver tolerances are
                         # used so that the results can be tested with reasonably tight
                         # tolerances without needing to match iteration counts (which
                         # could be affected by rounding errors, etc.).
                         "timestepping" => OptionsDict("type" => "KennedyCarpenterARK437",
                                                       "kinetic_ion_solver" => "full_explicit_ion_advance",
                                                       "nstep" => 10000,
                                                       "dt" => 7.071067811865475e-7,
                                                       "minimum_dt" => 7.071067811865474e-8,
                                                       "rtol" => 1.0e-11,
                                                       "atol" => 1.0e-14,
                                                       "nwrite" => 10000,
                                                       "exact_output_times" => true,
                                                       "high_precision_error_sum" => true),
                         "nonlinear_solver" => OptionsDict("nonlinear_max_iterations" => 100),
                         "r" => OptionsDict("ngrid" => 1,
                                            "nelement" => 1),
                         "z" => OptionsDict("ngrid" => 17,
                                            "nelement" => 16,
                                            "bc" => "periodic",
                                            "discretization" => "chebyshev_pseudospectral"),
                         "vpa" => OptionsDict("ngrid" => 6,
                                              "nelement" => 31,
                                              "L" => 20.784609690826528,
                                              "bc" => "zero",
                                              "discretization" => "chebyshev_pseudospectral"),
                         "vz" => OptionsDict("ngrid" => 6,
                                             "nelement" => 31,
                                             "L" => 20.784609690826528,
                                             "bc" => "zero",
                                             "discretization" => "chebyshev_pseudospectral"),
                         "ion_numerical_dissipation" => OptionsDict("force_minimum_pdf_value" => 0.0,
                                                                    "vpa_dissipation_coefficient" => 4.242640687119286),
                         "neutral_numerical_dissipation" => OptionsDict("force_minimum_pdf_value" => 0.0,
                                                                        "vz_dissipation_coefficient" => 0.42426406871192857),
                        )

if global_size[] > 2 && global_size[] % 2 == 0
    # Test using distributed-memory
    test_input["z"]["nelement_local"] = test_input["z"]["nelement"] ÷ 2
end

"""
Run a test for a single set of parameters
"""
function run_test(test_input, expected_p, expected_q, expected_vt; rtol=1.e-6,
                  qpar_atol=2.e-6, args...)
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
    p = undef
    q = undef
    vt = undef
    quietoutput() do
        # run simulation
        run_moment_kinetics(input)
    end

    if global_rank[] == 0
        quietoutput() do
            # Load and analyse output
            #########################

            path = joinpath(realpath(input["output"]["base_directory"]), name)

            # open the output file
            run_info = get_run_info_no_setup(path)

            pressure_zrt = get_variable(run_info, "electron_pressure")
            parallel_heat_flux_zrt = get_variable(run_info, "electron_parallel_heat_flux")
            thermal_speed_zrt = get_variable(run_info, "electron_thermal_speed")

            close_run_info(run_info)

            p = pressure_zrt[:,1,:]
            q = parallel_heat_flux_zrt[:,1,:]
            vt = thermal_speed_zrt[:,1,:]
        end

        # Regression test
        actual_p = p[begin:3:end, end]
        actual_q = q[begin:3:end, end]
        actual_vt = vt[begin:3:end, end]
        if expected_p == nothing
            # Error: no expected input provided
            println("data tested would be: ", actual_p)
            @test false
        else
            @test elementwise_isapprox(actual_p, expected_p, rtol=rtol, atol=0.0)
        end
        if expected_q == nothing
            # Error: no expected input provided
            println("data tested would be: ", actual_q)
            @test false
        else
            @test elementwise_isapprox(actual_q, expected_q, rtol=0.0, atol=qpar_atol)
        end
        if expected_vt == nothing
            # Error: no expected input provided
            println("data tested would be: ", actual_vt)
            @test false
        else
            @test elementwise_isapprox(actual_vt, expected_vt, rtol=rtol, atol=0.0)
        end
    end
end

function runtests()
    # Create a temporary directory for test output
    test_output_directory = get_MPI_tempdir()

    expected_p = [0.2912266585417873, 0.2902543624330032, 0.28792430978344075,
                  0.28548229473147724, 0.28384246425717863, 0.28316655242863065,
                  0.2829636736871053, 0.28233672371575724, 0.2818560684804787,
                  0.2819217850258229, 0.2822583755093484, 0.2823827721615315,
                  0.2828236538690061, 0.28397534851435324, 0.28576927524797635,
                  0.2875095272589636, 0.28823865551817063, 0.2890080800214039,
                  0.2912485879454332, 0.2945266724236888, 0.29778130426434,
                  0.2996483437262065, 0.3003187648990182, 0.30300274166441477,
                  0.30731924145835365, 0.3119315298636365, 0.31513680453994697,
                  0.3159841931420626, 0.31845475073366947, 0.3232154863021024,
                  0.3287724436519861, 0.333195499543554, 0.3348752062416145,
                  0.3365631172969006, 0.34108755833514826, 0.3469402297177915,
                  0.3521070366146737, 0.3548444063692133, 0.35579212239870905,
                  0.35941725850708117, 0.3647429903903947, 0.3698359230194339,
                  0.3730487420751043, 0.3738558988143583, 0.3761116794720789,
                  0.38006166546717024, 0.38403288142302383, 0.3867105490778209,
                  0.38761542294090523, 0.3884623070299775, 0.3904205952031905,
                  0.39226337604818395, 0.39321255665982724, 0.39344338708425086,
                  0.39347744373290805, 0.3933805180539749, 0.39253777459353145,
                  0.39085139006819514, 0.38927270817509996, 0.388805213880556,
                  0.38733115761835446, 0.3840629930279169, 0.37962417101340096,
                  0.37567650154540644, 0.37409141242592403, 0.3724543395463096,
                  0.36786264451518447, 0.36152744440454265, 0.3556083365746314,
                  0.35236008135635943, 0.351218280198776, 0.3467715470018111,
                  0.34001539481377363, 0.3332995336581819, 0.32892091438089555,
                  0.3278014020380251, 0.3246260079082095, 0.31886780380643875,
                  0.3127209531872688, 0.308249940529675, 0.3066465692484883,
                  0.3050866843362168, 0.3011542407539123, 0.29659122885303213,
                  0.29304173871331995, 0.29134086892961464]
    expected_q = [0.2355679155208151, 0.22887945534163132, 0.21041096151494268,
                  0.1854924635531075, 0.16262215758150048, 0.15017905091246722,
                  0.145816654170262, 0.1288565930423005, 0.10305658751455586,
                  0.07714920462486159, 0.05998325149536784, 0.05554611956523693,
                  0.04282964310320282, 0.0191517917951949, -0.007331771675645391,
                  -0.02769098913568107, -0.03528243080033772, -0.04284235161484274,
                  -0.06281020341147067, -0.08813221687009858, -0.11016352285791825,
                  -0.12176322110225359, -0.12577280987206277, -0.14110132103636192,
                  -0.16368642454986837, -0.1855217010639025, -0.19951542508566844,
                  -0.20306802474216107, -0.21309448275579607, -0.23109979670943276,
                  -0.2500569080988504, -0.2636173879505655, -0.2684140880931383,
                  -0.273036842635479, -0.2844340116829458, -0.29693104743529847,
                  -0.3056767476459368, -0.3093537993503738, -0.31046069578225644,
                  -0.31384845828793845, -0.3161072143169339, -0.3146312100975298,
                  -0.31143718003323645, -0.31030985313660653, -0.3063680812742302,
                  -0.29610649545791234, -0.27986766802342355, -0.26378295522195033,
                  -0.25695830623810584, -0.24970844070804382, -0.22835964589899121,
                  -0.19666308349901243, -0.1649036737549193, -0.14664890338942466,
                  -0.140097978768677, -0.11393871863004715, -0.07229712674584766,
                  -0.028765648823402816, 0.0006872908424582695, 0.008345664844609968,
                  0.03033692087257797, 0.07115370197051515, 0.11585277905040497,
                  0.14889622890500348, 0.16081299388669126, 0.17242312295801573,
                  0.20167354119441652, 0.23526702220291182, 0.26071798999931384,
                  0.27252902709915583, 0.2763401360833155, 0.28956422656946124,
                  0.3049448672748226, 0.31483731213392097, 0.3184381764490717,
                  0.3189966598786677, 0.3197691022597711, 0.3180338165404023,
                  0.3114571131880245, 0.3033021115996505, 0.29961405367073,
                  0.29560737451037744, 0.2834967271412086, 0.26518369476344644,
                  0.24682078800524532, 0.23632100055305744]
    expected_vt = [48.699287506193244, 48.66061223148292, 48.562473591185494,
                   48.44738750585768, 48.356893006442434, 48.31320367464671,
                   48.29877242828651, 48.246917828253615, 48.18065657782498,
                   48.129071741999454, 48.10304933992393, 48.097377745599154,
                   48.08352467822007, 48.0672376351584, 48.06375502762636,
                   48.07173585185883, 48.077098577116885, 48.08373343591118,
                   48.107500238010665, 48.15077649096474, 48.20055073283171,
                   48.2313699505044, 48.24277415080707, 48.289986829672884,
                   48.37029668539037, 48.46082660782321, 48.526058958176854,
                   48.54358265561334, 48.59528010133083, 48.697215775603354,
                   48.81957873384047, 48.91927429698015, 48.95763111014359,
                   48.996440701256795, 49.10175121379928, 49.24070705079878,
                   49.36598716579192, 49.43340637933757, 49.45692608193292,
                   49.54778998002585, 49.684142919437214, 49.81837423924499,
                   49.905475543045306, 49.927711141837264, 49.99072524427347,
                   50.1048379078491, 50.22642470702952, 50.314570235387016,
                   50.346068093277054, 50.37662669283179, 50.45311545885471,
                   50.54015858964053, 50.6054276734609, 50.63538607392922,
                   50.64498064589603, 50.677831644270555, 50.714216244632986,
                   50.73401504880432, 50.73767250707508, 50.73736741913412,
                   50.7336058469399, 50.715002851113404, 50.67589162985291,
                   50.63240476839586, 50.61316882225029, 50.592374082739276,
                   50.529567672673004, 50.43348611340068, 50.335038849542855,
                   50.2777068617995, 50.25701384418651, 50.17378452545543,
                   50.0393288891218, 49.895967819603705, 49.79708279454463,
                   49.77108639562322, 49.69571108156211, 49.5525231566822,
                   49.38951067947619, 49.26344624130706, 49.216506998753246,
                   49.169883043406216, 49.04771806575502, 48.89611621698131,
                   48.76869792938708, 48.70375705109887]

    @testset "Braginskii electron IMEX timestepping" verbose=use_verbose begin
        println("Braginskii electron IMEX timestepping tests")

        @testset "Split 3" begin
            test_input["output"]["base_directory"] = test_output_directory
            run_test(test_input, expected_p, expected_q, expected_vt)
        end
        @long @testset "Check other timestep - $type" for
                type ∈ ("KennedyCarpenterARK324",)

            timestep_check_input = deepcopy(test_input)
            timestep_check_input["output"]["base_directory"] = test_output_directory
            timestep_check_input["output"]["run_name"] = type
            timestep_check_input["timestepping"]["type"] = type
            run_test(timestep_check_input, expected_p, expected_q, expected_vt)
        end
    end

    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end # BraginskiiElectronsIMEX


using .BraginskiiElectronsIMEX

BraginskiiElectronsIMEX.runtests()
