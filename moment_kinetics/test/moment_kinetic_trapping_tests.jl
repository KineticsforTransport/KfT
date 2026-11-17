module moment_kinetic_trapping_tests

# Test generated from TOML input files

include("setup.jl")

using Base.Filesystem: tempname
using MPI

using moment_kinetics.interpolation: interpolate_to_grid_z
using moment_kinetics.load_data: get_run_info_no_setup, close_run_info,
                                 postproc_load_variable
using moment_kinetics.utils: merge_dict_with_kwargs!

# default inputs for tests
test_input = OptionsDict(
 "output" => OptionsDict(
     "display_timing_info" => true,
     "parallel_io" => true
    ),
 "r" => OptionsDict(
     "ngrid" => 1,
     "nelement" => 1
    ),
 "vperp" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 5,
     "L" => 10.0,
     "bc" => "zero"
    ),
 "evolve_moments" => OptionsDict(
     "pressure" => false,
     "density" => false,
     "parallel_flow" => false
    ),
 "ion_species_1" => OptionsDict(
     "initial_temperature" => 0.5,
     "initial_density" => 1.0
    ),
 "krook_collisions" => OptionsDict(
     "use_krook" => true,
     "frequency_option" => "reference_parameters"
    ),
 "vpa" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 20,
     "L" => 20.0,
     "element_spacing_option" => "coarse_tails",
     "bc" => "zero"
    ),
 "geometry" => OptionsDict(
     "pitch" => 1.0,
     "option" => "1D-mirror-STEP-edge-precise",
     "rhostar" => 1.0
    ),
 "z" => OptionsDict(
     "ngrid" => 5,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 32,
     "L" => 13.0,
     "element_spacing_option" => "sqrt",
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
     "source_T" => 0.1,
     "active" => true,
     "source_strength" => 1.0,
     "z_profile" => "wall_exp_decay",
     "z_width" => 0.5
    ),
 "ion_numerical_dissipation" => OptionsDict(
     "z_dissipation_coefficient" => 0.0
    ),
 "z_IC_ion_species_1" => OptionsDict(
     "initialization_option" => "gaussian",
     "density_amplitude" => 0.001,
     "temperature_amplitude" => 0.0,
     "density_phase" => 0.0,
     "upar_amplitude" => 1.0,
     "temperature_phase" => 0.0,
     "upar_phase" => 0.0
    ),
 "ion_source_1" => OptionsDict(
     "source_T" => 3.0,
     "active" => true,
     "source_strength" => 1.0,
     "z_profile" => "super_gaussian_4",
     "z_width" => 3.0
    ),
 "timestepping" => OptionsDict(
     "nstep" => 100,
     "steady_state_residual" => true,
     "dt" => 0.0005,
     "nwrite" => 1000,
     "type" => "SSPRK4",
     "nwrite_dfns" => 1000,
     "print_nT_live" => true
    )
)

test_input_2 = recursive_merge(test_input,
                               OptionsDict(
 "evolve_moments" => OptionsDict(
     "density" => true,
    )
))
test_input_3 = recursive_merge(test_input,
                               OptionsDict(
 "evolve_moments" => OptionsDict(
     "density" => true,
     "parallel_flow" => true
    )
))
test_input_4 = recursive_merge(test_input,
                               OptionsDict(
 "vpa" => OptionsDict(
     "nelement" => 40
    ),
 "timestepping" => OptionsDict(
     "nstep" => 200,
     "dt" => 0.00025,
     "nwrite" => 2000,
     "nwrite_dfns" => 2000
    ),
 "evolve_moments" => OptionsDict(
     "pressure" => true,
     "density" => true,
     "parallel_flow" => true
    )
))
# Here choose the names for each test
test_input = recursive_merge(test_input,
                               OptionsDict("output" => OptionsDict("run_name" => "nnn")))
test_input_2 = recursive_merge(test_input_2,
                               OptionsDict("output" => OptionsDict("run_name" => "ynn")))
test_input_3 = recursive_merge(test_input_3,
                               OptionsDict("output" => OptionsDict("run_name" => "yyn")))
test_input_4 = recursive_merge(test_input_4,
                               OptionsDict("output" => OptionsDict("run_name" => "yyy")))

"""
Run a test for a single set of parameters
"""
function run_test(test_input, expected_phi; rtol=0.0, atol=4.e-11, args...)
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

    @testset "moment_kinetic_trapping_tests tests" verbose=use_verbose begin
        println("moment_kinetic_trapping_tests tests")
        @testset "nnn" begin
            test_input["output"]["base_directory"] = test_output_directory
            run_test(test_input,
                     [-0.19324372338549253, -0.12394664909048042, -0.11848310763940267,
                     -0.11237655756044013, -0.10060528272993946, -0.085457852757343425,
                     -0.072325688861340234, -0.058482871340616006, -0.03830172381215019,
                     -0.015834323989597967, 0.0022401519554978858, 0.02061983540569344,
                     0.04603898095327362, 0.073080623961648244, 0.094606857668290981,
                     0.11671233623124865, 0.14644367298119826, 0.17370848737579603,
                     0.18995672320835286, 0.20121818354059876, 0.20927821616449582,
                     0.2117351028873686, 0.21123686634366276, 0.20839763443900564,
                     0.19899986870075595, 0.18158961752838651, 0.16327119975228566,
                     0.14157053205233042, 0.11138667940799137, 0.08262036973046856,
                     0.062015092280510387, 0.041773104715536941, 0.016112191185893011,
                     -0.0078537366090331677, -0.024977491111912789, -0.041681371898080771,
                     -0.062071901131539257, -0.079805684140684313, -0.091586412914311022,
                     -0.10255608959348979, -0.11447716212415673, -0.12100599771744516,
                     -0.13160220593576263])
        end
        @testset "ynn" begin
            test_input_2["output"]["base_directory"] = test_output_directory
            run_test(test_input_2,
                     [-0.19057601411404421, -0.12589157849774155, -0.11988668840291095,
                     -0.11380388386216991, -0.10168893643724496, -0.086292943132906058,
                     -0.072972831020543708, -0.058984670564465908, -0.038655889412640275,
                     -0.016086844873552567, 0.002039270108313737, 0.020455580683694189,
                     0.045912202044192094, 0.072987607270241547, 0.094539343094685535,
                     0.11665609078711901, 0.14639189455393042, 0.17366121280913521,
                     0.18990798512771523, 0.20116847096671575, 0.20922836574932255,
                     0.21169722802629798, 0.21119222830246509, 0.20834734459946314,
                     0.1989503625277072, 0.18154198173860653, 0.16322312017670451,
                     0.14151721402473516, 0.11133079609946517, 0.082538677229488489,
                     0.061908580940240945, 0.041640644797188257, 0.015939966843088585,
                     -0.0080809683427403734, -0.025265189053084836, -0.042055809862658063,
                     -0.062607380158838863, -0.080552121840630403, -0.092509998576281519,
                     -0.10373812517375017, -0.11579179384670518, -0.12279525858473683,
                     -0.13292592376654985])
        end
        @testset "yyn" begin
            test_input_3["output"]["base_directory"] = test_output_directory
            run_test(test_input_3,
                     [-0.1854411880917714, -0.12544059684898487, -0.12010759862631813,
                      -0.11401348790768276, -0.10175770324861458, -0.086305872014501,
                      -0.07297851729577486, -0.05898009463640507, -0.038651508338221445,
                      -0.016098918189006124, 0.0020013379882024124, 0.02051370477667807,
                      0.04590487080104194, 0.07297542044186256, 0.09455289881447843,
                      0.11666593733534883, 0.14638803850742918, 0.17358703657410568,
                      0.189907873847179, 0.20116257917301048, 0.2092252647041966,
                      0.21168408654495938, 0.21118323749796497, 0.20831512723267231,
                      0.19897300236277612, 0.18164305533525732, 0.16319947848656977,
                      0.1415153380291998, 0.11136058199363537, 0.08250741178535237,
                      0.06191526924758056, 0.041636867646829376, 0.01597433219687073,
                      -0.008099984754525419, -0.02526863752380065, -0.042050289462751356,
                      -0.06260479427755378, -0.08056438173339854, -0.0925608935252589,
                      -0.10371144727825357, -0.11615631729263534, -0.12284151294565937,
                      -0.13237086114783292])
        end
        @testset "yyy" begin
            test_input_4["output"]["base_directory"] = test_output_directory
            run_test(test_input_4,
                     [-0.19276271645096404, -0.12665828933624915, -0.12045506239788012,
                      -0.11428757028711951, -0.10221870621888118, -0.08680924587186627,
                      -0.07352771872842631, -0.05959041994335114, -0.039316349893317046,
                      -0.016558585365107438, 0.001533751704333245, 0.01997903691419602,
                      0.04542597161764525, 0.07247775384869874, 0.0940259844143765,
                      0.1161255124879753, 0.14585483563857468, 0.17312934833909033,
                      0.18938853724749818, 0.20065999779250254, 0.20872792442677582,
                      0.21118795777766813, 0.21068792430986127, 0.2078436828819743,
                      0.19852503393517015, 0.18101226845262972, 0.1626852829693904,
                      0.1409795845263093, 0.11080138727508962, 0.08203061816808255,
                      0.061408831867401606, 0.04115571070863462, 0.015388458279860224,
                      -0.008505293350487752, -0.025496400022871873, -0.04273762860364659,
                      -0.06320241247021074, -0.0810843002257014, -0.09305980747231477,
                      -0.1042169689612063, -0.11633412011279172, -0.12314134254967682,
                      -0.13282179793154023])
        end
    end
    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end

using .moment_kinetic_trapping_tests

moment_kinetic_trapping_tests.runtests()
