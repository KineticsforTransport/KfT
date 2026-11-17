module SoundWaveTests

include("setup.jl")

using Base.Filesystem: tempname
#using Plots: plot, plot!, gui

using moment_kinetics.array_allocation: allocate_float
using moment_kinetics.input_structs: netcdf
using moment_kinetics.file_io: io_has_implementation
using moment_kinetics.load_data: open_readonly_output_file
using moment_kinetics.load_data: load_fields_data, load_time_data
using moment_kinetics.load_data: load_species_data, load_coordinate_data
using moment_kinetics.analysis: analyze_fields_data
using moment_kinetics.analysis: fit_delta_phi_mode
using moment_kinetics.utils: merge_dict_with_kwargs!

const analytical_rtol = 3.e-2
const regression_rtol = 1.e-10
const regression_range = 5:10

# Use "netcdf" to test the NetCDF I/O if it is available (or if we are forcing optional
# dependencies to be used, e.g. for CI tests), otherwise fall back to "hdf5".
const binary_format = (force_optional_dependencies || io_has_implementation(Val(netcdf))) ?
                      "netcdf" : "hdf5"

# default inputs for tests
test_input_finite_difference = OptionsDict("composition" => OptionsDict("n_ion_species" => 1,
                                                                        "n_neutral_species" => 1,
                                                                        "electron_physics" => "boltzmann_electron_response",
                                                                        "T_e" => 1.0),
                                           "ion_species_1" => OptionsDict("initial_density" => 0.5,
                                                                          "initial_temperature" => 0.3333333333333333),
                                           "z_IC_ion_species_1" => OptionsDict("initialization_option" => "sinusoid",
                                                                               "density_amplitude" => 0.001,
                                                                               "density_phase" => 0.0,
                                                                               "upar_amplitude" => 0.0,
                                                                               "upar_phase" => 0.0,
                                                                               "temperature_amplitude" => 0.0,
                                                                               "temperature_phase" => 0.0),
                                           "neutral_species_1" => OptionsDict("initial_density" => 0.5,
                                                                              "initial_temperature" => 0.3333333333333333),
                                           "z_IC_neutral_species_1" => OptionsDict("initialization_option" => "sinusoid",
                                                                                   "density_amplitude" => 0.001,
                                                                                   "density_phase" => 0.0,
                                                                                   "upar_amplitude" => 0.0,
                                                                                   "upar_phase" => 0.0,
                                                                                   "temperature_amplitude" => 0.0,
                                                                                   "temperature_phase" => 0.0),
                                           "output" => OptionsDict("run_name" => "finite_difference",
                                                                   "binary_format" => binary_format),
                                           "evolve_moments" => OptionsDict("density" => false,
                                                                           "parallel_flow" => false,
                                                                           "pressure" => false),
                                           "reactions" => OptionsDict("charge_exchange_frequency" => 0.8885765876316732,
                                                                      "ionization_frequency" => 0.0),
                                           "timestepping" => OptionsDict("nstep" => 1500,
                                                                         "dt" => 0.001414213562373095,
                                                                         "nwrite" => 20,
                                                                         "split_operators" => false),
                                           "r" => OptionsDict("ngrid" => 1,
                                                              "nelement" => 1,
                                                              "discretization" => "finite_difference"),
                                           "z" => OptionsDict("ngrid" => 100,
                                                              "nelement" => 1,
                                                              "bc" => "periodic",
                                                              "discretization" => "finite_difference"),
                                           "vperp" => OptionsDict("ngrid" => 1,
                                                                  "nelement" => 1,
                                                                  "L" => 1.4142135623730951,
                                                                  "discretization" => "finite_difference"),
                                           "vpa" => OptionsDict("ngrid" => 180,
                                                                "nelement" => 1,
                                                                "L" => 11.313708498984761,
                                                                "bc" => "periodic",
                                                                "discretization" => "finite_difference"),
                                           "vz" => OptionsDict("ngrid" => 180,
                                                               "nelement" => 1,
                                                               "L" => 11.313708498984761,
                                                               "bc" => "periodic",
                                                               "discretization" => "finite_difference"),
                                          )

test_input_finite_difference_split_1_moment =
    recursive_merge(test_input_finite_difference,
                    OptionsDict("output" => OptionsDict("run_name" => "finite_difference_split_1_moment"),
                                "evolve_moments" => OptionsDict("density" => true))
                   )

test_input_finite_difference_split_2_moments =
    recursive_merge(test_input_finite_difference_split_1_moment,
                    OptionsDict("output" => OptionsDict("run_name" => "finite_difference_split_2_moments"),
                                "evolve_moments" => OptionsDict("parallel_flow" => true),
                                "vpa" => OptionsDict("ngrid" => 270, "L" => 16.970562748477143),
                                "vz" => OptionsDict("ngrid" => 270, "L" => 16.970562748477143))
                   )

test_input_finite_difference_split_3_moments =
    recursive_merge(test_input_finite_difference_split_2_moments,
                    OptionsDict("output" => OptionsDict("run_name" => "finite_difference_split_3_moments"),
                                "evolve_moments" => OptionsDict("pressure" => true),
                                "vpa" => OptionsDict("ngrid" => 270, "L" => 20.784609690826528),
                                "vz" => OptionsDict("ngrid" => 270, "L" => 20.784609690826528))
                   )

test_input_chebyshev = recursive_merge(test_input_finite_difference,
                                       OptionsDict("output" => OptionsDict("run_name" => "chebyshev_pseudospectral"),
                                                   "z" => OptionsDict("discretization" => "chebyshev_pseudospectral",
                                                                      "ngrid" => 9,
                                                                      "nelement" => 2),
                                                   "vpa" => OptionsDict("discretization" => "chebyshev_pseudospectral",
                                                                        "ngrid" => 17,
                                                                        "nelement" => 8),
                                                   "vz" => OptionsDict("discretization" => "chebyshev_pseudospectral",
                                                                       "ngrid" => 17,
                                                                       "nelement" => 8),
                                                  ))

test_input_chebyshev_split_1_moment =
    recursive_merge(test_input_chebyshev,
                    OptionsDict("output" => OptionsDict("run_name" => "chebyshev_pseudospectral_split_1_moment"),
                                "evolve_moments" => OptionsDict("density" => true)))

test_input_chebyshev_split_2_moments =
    recursive_merge(test_input_chebyshev_split_1_moment,
                    OptionsDict("output" => OptionsDict("run_name" => "chebyshev_pseudospectral_split_2_moments"),
                                "evolve_moments" => OptionsDict("parallel_flow" => true)))

test_input_chebyshev_split_3_moments =
    recursive_merge(test_input_chebyshev_split_2_moments,
                    OptionsDict("output" => OptionsDict("run_name" => "chebyshev_pseudospectral_split_3_moments"),
                                "evolve_moments" => OptionsDict("pressure" => true),
                                "vpa" => OptionsDict("L" => 13.856406460551018),
                                "vz" => OptionsDict("L" => 13.856406460551018)))


"""
Run a sound-wave test for a single set of parameters
"""
# Note 'name' should not be shared by any two tests in this file
function run_test(test_input, analytic_frequency, analytic_growth_rate,
                  regression_phi, itime_min=50; args...)
    # by passing keyword arguments to run_test, args becomes a Dict which can be used to
    # update the default inputs

    # Make a copy to make sure nothing modifies the input Dicts defined in this test
    # script.
    input = deepcopy(test_input)

    # Convert keyword arguments to a unique name
    function stringify_arg(key, value)
        if isa(value, AbstractDict)
            return string(string(key)[1], (stringify_arg(k, v) for (k, v) in value)...)
        else
            return string(string(key)[1], value)
        end
    end
    name = input["output"]["run_name"]
    shortname = name
    if length(args) > 0
        name = string(name, "_", (string(k, "-", v, "_") for (k, v) in args)...)
        shortname = string(shortname, "_", (stringify_arg(k, v) for (k, v) in args)...)

        # Remove trailing "_"
        name = chop(name)
    end

    # Provide some progress info
    println("    - testing ", name)

    # Update default inputs with values to be changed
    merge_dict_with_kwargs!(input; args...)
    input["output"]["run_name"] = shortname

    # Suppress console output while running
    phi_fit = undef
    phi = undef
    quietoutput() do
        # run simulation
        run_moment_kinetics(input)
    end

    if global_rank[] == 0
        quietoutput() do

            # Load and analyse output
            #########################

            path = joinpath(realpath(input["output"]["base_directory"]), shortname, shortname)

            # open the netcdf file and give it the handle 'fid'
            fid = open_readonly_output_file(path,"moments")

            # load space-time coordinate data
            z, z_spectral = load_coordinate_data(fid, "z"; ignore_MPI=true)
            r, r_spectral = load_coordinate_data(fid, "r"; ignore_MPI=true)
            n_ion_species, n_neutral_species = load_species_data(fid)
            ntime, time = load_time_data(fid)
            
            # load fields data
            phi_zrt, Er_zrt, Ez_zrt = load_fields_data(fid)

            close(fid)
            
            ir0 = 1 
            
            phi = phi_zrt[:,ir0,:]

            # analyze the fields data
            phi_fldline_avg, delta_phi = analyze_fields_data(phi, ntime, z)

            # use a fit to calculate the damping rate and growth rate of the perturbed
            # electrostatic potential
            itime_max = ntime
            iz0 = cld(z.n, 3)
            shifted_time = allocate_float(ntime)
            @. shifted_time = time - time[itime_min]
            @views phi_fit = fit_delta_phi_mode(shifted_time[itime_min:itime_max], z.grid,
                                                delta_phi[:, itime_min:itime_max])
            ## The following plot code (copied from post_processing.jl) may be helpful for
            ## debugging tests. Uncomment to use, and also uncomment
            ## `using Plots: plot, plot!, gui at the top of the file.
            #L = z.grid[end] - z.grid[begin]
            #fitted_delta_phi =
            #    @. (phi_fit.amplitude0 * cos(2.0 * π * (z.grid[iz0] + phi_fit.offset0) / L)
            #        * exp(phi_fit.growth_rate * shifted_time)
            #        * cos(phi_fit.frequency * shifted_time + phi_fit.phase))
            #@views plot(time, abs.(delta_phi[iz0,:]), xlabel="t*z.L/vti", ylabel="δϕ", yaxis=:log)
            #plot!(time, abs.(fitted_delta_phi))
            #gui()
        end

        # Check the fit errors are not too large, otherwise we are testing junk
        @test phi_fit.amplitude_fit_error < 1.e-1
        @test phi_fit.offset_fit_error < 5.e-6
        @test phi_fit.cosine_fit_error < 5.e-8

        # analytic_frequency and analytic_growth rate are the analytically expected values
        # (from F. Parra's calculation).
        @test isapprox(phi_fit.frequency, analytic_frequency, rtol=analytical_rtol)
        @test isapprox(phi_fit.growth_rate, analytic_growth_rate, rtol=analytical_rtol)

        # Test some values of phi for a regression test, which can use with tighter
        # tolerances than the analytic test.
        if regression_phi === nothing
            println("Values tested would be: ", phi[regression_range,end])
            @test false
        else
            @test isapprox(phi[regression_range,end], regression_phi, rtol=regression_rtol)
        end
    end
end


# run_test_set_* functions call run_test for various parameters, and record the expected
# values to be used for regression_frequency and regression_growth_rate for each
# particular case

function run_test_set_finite_difference()
    #n_i=n_n, T_e=1
    @long run_test(test_input_finite_difference, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.6931471878909734, -0.693147187702995, -0.6931471875145545,
                    -0.6931471873358439, -0.6931471871575013, -0.6931471869618118];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_finite_difference, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471896500133, -0.6931471895025301, -0.6931471893471256,
              -0.6931471891899375, -0.6931471890256098, -0.6931471888439447])
    @long run_test(test_input_finite_difference, 2*π*0.0 * sqrt(2), -2*π*0.3235 * sqrt(2),
                   [-0.6931469705307232, -0.6931469744636533, -0.6931469792610312,
                    -0.6931469849086842, -0.6931469913785957, -0.6931469986348702];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    @long run_test(test_input_finite_difference, 2*π*0.0 * sqrt(2), -2*π*0.2963 * sqrt(2),
                   [-0.6931468579689712, -0.6931468639804542, -0.6931468713014981,
                    -0.6931468799077714, -0.6931468897592976, -0.6931469008064994];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    @long run_test(test_input_finite_difference, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010001233181730678, -0.00010001214384163556,
                    -0.00010001195540095647, -0.00010001177668924657,
                    -0.00010001159834028284, -0.00010001140264423694];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    @long run_test(test_input_finite_difference, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010001221254433313, -0.00010001202677078748,
                    -0.00010001184100445905, -0.0001000116654257788,
                    -0.00010001149066097346, -0.00010001129898267139];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    @long run_test(test_input_finite_difference, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
                   [-9.210340384403116, -9.210340384323528, -9.210340384234046,
                    -9.210340384136575, -9.21034038402964, -9.21034038391042];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    @long run_test(test_input_finite_difference, 2*π*0.0 * sqrt(2), -2*π*0.5112 * sqrt(2),
                   [-9.210340471125871, -9.210340472461294, -9.210340474082226,
                    -9.210340475982138, -9.210340478153377, -9.210340480587199];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    @long run_test(test_input_finite_difference, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
                   [-0.34657359201449195, -0.34657359197740367, -0.34657359194126197,
                    -0.346573591903444, -0.34657359185860664, -0.34657359180391195], 30;
                   composition = OptionsDict("T_e" => 0.5), 
                   timestepping = OptionsDict("nstep" => 1300),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    @long run_test(test_input_finite_difference, 2*π*0.0 * sqrt(2), -2*π*0.2727 * sqrt(2),
                   [-0.3465733721093574, -0.3465733761667264, -0.3465733811031762,
                    -0.3465733869008616, -0.3465733935343765, -0.3465734009730069];
                   composition = OptionsDict("T_e" => 0.5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    @long run_test(test_input_finite_difference, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
                   [-2.772608882427265, -2.7726085100293045, -2.772608058199096,
                    -2.7726075287624155, -2.772606923796338, -2.772606245649231];
                   composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_finite_difference_split_1_moment()
    #n_i=n_n, T_e=1
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
             [-0.6931471878126522, -0.693147187632401, -0.6931471874480738,
              -0.6931471872679154, -0.6931471870864971, -0.6931471868905795];
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471867765525, -0.6931471866343989, -0.6931471864823142,
              -0.6931471863254078, -0.693147186160579, -0.6931471859802977])
    run_test(test_input_finite_difference_split_1_moment, 2*π*0.0, -2*π*0.3235 * sqrt(2),
             [-0.6931469658591621, -0.6931469697975196, -0.6931469745986458,
              -0.6931469802474806, -0.6931469867180096, -0.6931469939766792];
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    run_test(test_input_finite_difference_split_1_moment, 2*π*0.0, -2*π*0.2963 * sqrt(2),
             [-0.6931468532790588, -0.6931468592963252, -0.6931468666215331,
              -0.6931468752294686, -0.6931468850821589, -0.6931468961323719];
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
             [-0.00010001225310683516, -0.00010001207285953319, -0.00010001188853164252,
              -0.00010001170837172392, -0.00010001152694746793, -0.00010001133102480287];
             ion_species_1 = OptionsDict("initial_density" => 0.9999), 
             neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
             [-0.00010001213374381343, -0.00010001195569730463, -0.00010001177404509702,
              -0.0001000115970189853, -0.0001000114191805532, -0.0001000112272744106];
             ion_species_1 = OptionsDict("initial_density" => 0.9999), 
             neutral_species_1 = OptionsDict("initial_density" => 0.0001),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
             [-9.210340375710334, -9.210340375633665, -9.210340375546723,
              -9.210340375451267, -9.210340375346378, -9.210340375230128];
             ion_species_1 = OptionsDict("initial_density" => 0.0001), 
             neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    run_test(test_input_finite_difference_split_1_moment, 2*π*0.0, -2*π*0.5112 * sqrt(2),
             [-9.210340299661922, -9.210340300998629, -9.21034030262106,
              -9.210340304522843, -9.210340306696159, -9.210340309132416];
             ion_species_1 = OptionsDict("initial_density" => 0.0001), 
             neutral_species_1 = OptionsDict("initial_density" => 0.9999),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
             [-0.34657359193147635, -0.34657359189694775, -0.34657359186157893,
              -0.3465735918241752, -0.346573591781001, -0.3465735917294017], 30;
             composition = OptionsDict("T_e" => 0.5),
             timestepping = OptionsDict("nstep" => 1300),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_finite_difference_split_1_moment, 2*π*0.0, -2*π*0.2727 * sqrt(2),
             [-0.34657336973345404, -0.3465733737937447, -0.3465733787327239,
              -0.34657338453223374, -0.34657339116759117, -0.3465733986089023];
             composition = OptionsDict("T_e" => 0.5),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    run_test(test_input_finite_difference_split_1_moment, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
             [-2.772608874435827, -2.772608502007103, -2.7726080501068613,
              -2.7726075205725347, -2.772606915491342, -2.772606237232529];
              composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_finite_difference_split_2_moments()
    #n_i=n_n, T_e=1
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
             [-0.6931471880976617, -0.6931471879224558, -0.6931471877390526,
              -0.6931471875549782, -0.6931471873623963, -0.6931471871467376];
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471870142223, -0.6931471868675833, -0.6931471867121105,
              -0.6931471865544573, -0.6931471863882103, -0.6931471862008788])
    run_test(test_input_finite_difference_split_2_moments, 2*π*0.0, -2*π*0.3235 * sqrt(2),
             [-0.6931469660834078, -0.6931469700165075, -0.6931469748143118,
              -0.6931469804632058, -0.6931469869339888, -0.6931469941867228];
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    run_test(test_input_finite_difference_split_2_moments, 2*π*0.0, -2*π*0.2963 * sqrt(2),
             [-0.6931468535072797, -0.6931468595196247, -0.6931468668415703,
              -0.6931468754493064, -0.6931468853017455, -0.6931468963452997];
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
             [-0.00010001253810727225, -0.00010001236290277446, -0.00010001217949812652,
              -0.00010001199542261477, -0.00010001180283828002, -0.00010001158717631795];
             ion_species_1 = OptionsDict("initial_density" => 0.9999), 
             neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
             [-0.00010001241878877488, -0.000100012245787513, -0.00010001206505788191,
              -0.00010001188411873082, -0.00010001169511922065, -0.00010001148347711205];
             ion_species_1 = OptionsDict("initial_density" => 0.9999), 
             neutral_species_1 = OptionsDict("initial_density" => 0.0001),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
             [-9.210340375902973, -9.210340375812587, -9.21034037571806,
              -9.210340375624837, -9.210340375527633, -9.210340375415985];
             ion_species_1 = OptionsDict("initial_density" => 0.0001), 
             neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    run_test(test_input_finite_difference_split_2_moments, 2*π*0.0, -2*π*0.5112 * sqrt(2),
             [-9.210340299844662, -9.210340301166584, -9.210340302781027,
              -9.210340304685795, -9.210340306868739, -9.210340309310785];
             ion_species_1 = OptionsDict("initial_density" => 0.0001), 
             neutral_species_1 = OptionsDict("initial_density" => 0.9999),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
             [-0.3465735917666568, -0.34657359174673213, -0.3465735917228738,
              -0.3465735916993949, -0.3465735916797001, -0.3465735916641262], 30;
             composition = OptionsDict("T_e" => 0.5),
             timestepping = OptionsDict("nstep" => 1300), z = OptionsDict("ngrid" => 150),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_finite_difference_split_2_moments, 2*π*0.0, -2*π*0.2727 * sqrt(2),
             [-0.3465733699053265, -0.34657337396442434, -0.34657337890186546,
              -0.34657338470094173, -0.34657339133462495, -0.3465733987683617];
             composition = OptionsDict("T_e" => 0.5),
             reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    run_test(test_input_finite_difference_split_2_moments, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
             [-2.7726088854800315, -2.7726085130196982, -2.7726080609237274,
              -2.7726075310187794, -2.7726069254154453, -2.772606246520561];
             composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_finite_difference_split_3_moments()
    #n_i=n_n, T_e=1
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.6931471880707439, -0.6931471879690402, -0.6931471878177301,
                    -0.6931471876182134, -0.6931471873673457, -0.6931471870645756];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_finite_difference_split_3_moments, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471871103377, -0.6931471870343805, -0.6931471869090101,
              -0.6931471867365716, -0.693147186514809, -0.6931471862424651])
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*0.0, -2*π*0.3235 * sqrt(2),
                   [-0.6931469662670939, -0.693146970267982, -0.6931469750959321,
                    -0.6931469807337137, -0.6931469871546665, -0.6931469943276766];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*0.0, -2*π*0.2963 * sqrt(2),
                   [-0.6931468536986971, -0.6931468597793667, -0.6931468671315708,
                    -0.6931468757278696, -0.693146885529709, -0.6931468964924642];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.0001000125112140929, -0.00010001240950838069,
                    -0.00010001225819871616, -0.00010001205867845754,
                    -0.00010001180781057802, -0.00010001150503635756];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010001239192912761, -0.0001000122924263182,
                    -0.00010001214379167052, -0.00010001194740788359,
                    -0.00010001170012394042, -0.00010001140136746379];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
                   [-9.21034037612422, -9.21034037610433, -9.210340376039078,
                    -9.210340375931095, -9.210340375778818, -9.210340375580639];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*0.0, -2*π*0.5112 * sqrt(2),
                   [-9.21034030017531, -9.210340301568875, -9.210340303212021,
                    -9.210340305100555, -9.210340307225575, -9.210340309576397];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
                   [-0.34657359187283704, -0.34657359183783815, -0.3465735918120044,
                    -0.3465735917907561, -0.3465735917644153, -0.3465735917250847], 30;
                   composition = OptionsDict("T_e" => 0.5),
                   timestepping = OptionsDict("nstep" => 1300),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*0.0, -2*π*0.2727 * sqrt(2),
                   [-0.34657337010004446, -0.34657337420743617, -0.34657337916529274,
                    -0.3465733849553589, -0.3465733915522979, -0.3465733989268312];
                   composition = OptionsDict("T_e" => 0.5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    @long run_test(test_input_finite_difference_split_3_moments, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
                   [-2.7726088930013004, -2.772608520437762, -2.7726080681841303,
                    -2.77260753805663, -2.7726069321723914, -2.7726062529596893],;
                   composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_chebyshev()
    #n_i=n_n, T_e=1
    @long run_test(test_input_chebyshev, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.6931471805961407, -0.6931471765630876, -0.6931471742786531,
                    -0.6931471736212292, -0.693147173557635, -0.6931471736212295];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_chebyshev, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471805805466, -0.69314717719617, -0.6931471752710278,
              -0.6931471747226033, -0.6931471746642256, -0.6931471747226029])
    @long run_test(test_input_chebyshev, 2*π*0.0, -2*π*0.3235 * sqrt(2),
                   [-0.6931471805855942, -0.6931473062881018, -0.6931473797751034,
                    -0.6931474013292136, -0.6931474028946445, -0.693147401329213];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    @long run_test(test_input_chebyshev, 2*π*0.0, -2*π*0.2963 * sqrt(2),
                   [-0.6931471805861977, -0.6931473720683567, -0.6931474839886861,
                    -0.6931475168064498, -0.6931475192016895, -0.6931475168064497];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    @long run_test(test_input_chebyshev, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000503652579628, -0.00010000100344075166, -9.999871898506605e-5,
                    -9.999806155560843e-5, -9.999799796067053e-5, -9.999806155571947e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    @long run_test(test_input_chebyshev, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000503652835005, -0.00010000107315673933, -9.999882942888377e-5,
                    -9.999818393662746e-5, -9.999812121907552e-5, -9.999818393518403e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    @long run_test(test_input_chebyshev, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
                   [-9.210340371980537, -9.210340370097294, -9.210340369029236,
                    -9.210340368735709, -9.210340368693467, -9.210340368735709];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    @long run_test(test_input_chebyshev, 2*π*0.0, -2*π*0.5112 * sqrt(2),
                   [-9.21034037197519, -9.210340414554372, -9.21034043945508,
                    -9.210340446771436, -9.210340447284608, -9.210340446771436];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    @long run_test(test_input_chebyshev, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
                   [-0.34657359028138823, -0.3465735895309023, -0.3465735890566611,
                    -0.3465735889775626, -0.34657358897367324, -0.34657358897756235],
                   30; composition = OptionsDict("T_e" => 0.5),
                   timestepping = OptionsDict("nstep" => 1300),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    @long run_test(test_input_chebyshev, 2*π*0.0, -2*π*0.2727 * sqrt(2),
                   [-0.34657359029097584, -0.34657371934209186, -0.34657379475484945,
                    -0.346573816865366, -0.3465738184803065, -0.3465738168653658];
                   composition = OptionsDict("T_e" => 0.5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    @long run_test(test_input_chebyshev, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
                   [-2.772588723097548, -2.7725769408799374, -2.772570056405588,
                    -2.7725680387868135, -2.772567890042227, -2.7725680387868143];
                   composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_chebyshev_split_1_moment()
    #n_i=n_n, T_e=1
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.6931471798887808, -0.693147177311798, -0.6931471734976317,
                    -0.693147174320916, -0.6931471728800184, -0.6931471743209164];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_chebyshev_split_1_moment, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471798714258, -0.6931471779264865, -0.6931471745185674,
              -0.6931471754126449, -0.6931471739811726, -0.6931471754126449])
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*0.0, -2*π*0.3235 * sqrt(2),
                   [-0.6931471798770498, -0.6931473070165636, -0.6931473790257755,
                    -0.6931474020177143, -0.6931474022169881, -0.6931474020177147];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*0.0, -2*π*0.2963 * sqrt(2),
                   [-0.693147179877793, -0.6931473727963409, -0.6931474832383637,
                    -0.6931475174937916, -0.6931475185232685, -0.6931475174937918];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000432916594435, -0.0001000017521459994, -9.999793796787836e-5,
                    -9.999876123882836e-5, -9.999732034123214e-5, -9.999876123927249e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000432916783192, -0.00010000182186143193, -9.999804841325047e-5,
                    -9.999888362095782e-5, -9.999744360196874e-5, -9.999888361984747e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
                   [-9.21034037126966, -9.210340370808568, -9.210340368306326,
                    -9.210340369415363, -9.210340368003612, -9.210340369415361];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*0.0, -2*π*0.5112 * sqrt(2),
                   [-9.210340371263273, -9.21034041525851, -9.210340438745584,
                    -9.210340447449697, -9.210340446598197, -9.210340447449697];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
                   [-0.34657358993342785, -0.34657358985387404, -0.34657358873350347,
                    -0.34657358929201915, -0.3465735885789893, -0.3465735892920189],
                   30; composition = OptionsDict("T_e" => 0.5),
                   timestepping = OptionsDict("nstep" => 1300),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*0.0, -2*π*0.2727 * sqrt(2),
                   [-0.3465735899333066, -0.34657371970143575, -0.34657379438589264,
                    -0.34657381720435393, -0.3465738181404627, -0.34657381720435393];
                   composition = OptionsDict("T_e" => 0.5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    @long run_test(test_input_chebyshev_split_1_moment, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
                   [-2.7725887202721675, -2.7725769438619237, -2.772570053870096,
                    -2.772568041774934, -2.7725678875212694, -2.772568041774933];
                   composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_chebyshev_split_2_moments()
    #n_i=n_n, T_e=1
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.6931471798501451, -0.6931471773010914, -0.6931471735192708,
                    -0.6931471743562659, -0.6931471729237529, -0.6931471743562659];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_chebyshev_split_2_moments, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471798333434, -0.6931471779184996, -0.6931471745406236,
              -0.6931471754430353, -0.6931471740189036, -0.6931471754430355])
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*0.0, -2*π*0.3235 * sqrt(2),
                   [-0.693147179841931, -0.6931473070098627, -0.6931473790459587,
                    -0.6931474020440795, -0.6931474022498966, -0.6931474020440795];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*0.0, -2*π*0.2963 * sqrt(2),
                   [-0.6931471798427714, -0.6931473727895734, -0.6931474832583618,
                    -0.6931475175199603, -0.6931475185559783, -0.6931475175199603];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000429052765156, -0.00010000174143738261,
                    -9.99979596062904e-5, -9.999879658331498e-5, -9.999736407318211e-5,
                    -9.999879658431427e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000429053720043, -0.00010000181115747854,
                    -9.999807005310595e-5, -9.999891896622167e-5, -9.999748733680559e-5,
                    -9.999891896622167e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
                   [-9.210340371231817, -9.210340370802884, -9.210340368328923,
                    -9.21034036944169, -9.210340368036372, -9.210340369441692];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*0.0, -2*π*0.5112 * sqrt(2),
                   [-9.210340371225158, -9.210340415253803, -9.210340438768354,
                    -9.210340447473778, -9.21034044662836, -9.210340447473778];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
                   [-0.34657358990573656, -0.3465735898502191, -0.34657358875039707,
                    -0.346573589308585, -0.3465735885985403, -0.34657358930858523],
                   40; composition = OptionsDict("T_e" => 0.5),
                   timestepping = OptionsDict("nstep" => 1300, "nwrite" => 10),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0 * sqrt(2)))
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*0.0, -2*π*0.2727 * sqrt(2),
                   [-0.3465735899090165, -0.34657371969765555, -0.3465737944001191,
                    -0.34657381722103725, -0.34657381816134814, -0.34657381722103753];
                   composition = OptionsDict("T_e" => 0.5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    @long run_test(test_input_chebyshev_split_2_moments, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
                   [-2.772588720085524, -2.7725769438106878, -2.7725700538858162,
                    -2.772568041738153, -2.772567887524583, -2.772568041738153];
                   composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function run_test_set_chebyshev_split_3_moments()
    #n_i=n_n, T_e=1
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.6931471799373886, -0.6931471773221646, -0.6931471735214639,
                    -0.6931471743234282, -0.693147172875601, -0.6931471743234282];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    run_test(test_input_chebyshev_split_3_moments, 2*π*1.4240 * sqrt(2), -2*π*0.6379 * sqrt(2),
             [-0.6931471799082319, -0.6931471779413488, -0.6931471745498324,
              -0.6931471753838919, -0.6931471739471132, -0.6931471753838928])
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*0.0, -2*π*0.3235 * sqrt(2),
                   [-0.6931471799063286, -0.6931473070400977, -0.693147379076884,
                    -0.6931474019869261, -0.6931474021820054, -0.6931474019869253];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*1.8 * sqrt(2)))
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*0.0, -2*π*0.2963 * sqrt(2),
                   [-0.6931471799070094, -0.693147372821742, -0.6931474832928677,
                    -0.6931475174663773, -0.6931475184916174, -0.6931475174663768];
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i>>n_n T_e=1
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000437776837497, -0.00010000176251707462,
                    -9.999796180253116e-5, -9.999876374796125e-5, -9.999731591732818e-5,
                    -9.999876374751711e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001))
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*1.4467 * sqrt(2), -2*π*0.6020 * sqrt(2),
                   [-0.00010000437777037358, -0.00010000183223272923,
                    -9.999807225101221e-5, -9.999888612420593e-5, -9.999743917562206e-5,
                    -9.999888612287352e-5];
                   ion_species_1 = OptionsDict("initial_density" => 0.9999), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.0001),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i<<n_n T_e=1
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*1.3954 * sqrt(2), -2*π*0.6815 * sqrt(2),
                   [-9.21034037129687, -9.21034037082684, -9.210340368343163,
                    -9.210340369359246, -9.210340367943102, -9.210340369359244];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999))
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*0.0, -2*π*0.5112 * sqrt(2),
                   [-9.2103403712863, -9.210340415282188, -9.210340438792905,
                    -9.210340447390209, -9.210340446535636, -9.210340447390207];
                   ion_species_1 = OptionsDict("initial_density" => 0.0001), 
                   neutral_species_1 = OptionsDict("initial_density" => 0.9999),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=0.5
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*1.2671 * sqrt(2), -2*π*0.8033 * sqrt(2),
                   [-0.34657358995311743, -0.34657358987978654, -0.3465735887122216,
                    -0.3465735892962168, -0.3465735885613302, -0.34657358929621646],
                   80; composition = OptionsDict("T_e" => 0.5),
                   timestepping = OptionsDict("nstep" => 1300, "nwrite" => 5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*0.0))
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*0.0, -2*π*0.2727 * sqrt(2),
                   [-0.3465735899482742, -0.3465737197177559, -0.34657379442476904,
                    -0.34657381717123564, -0.34657381810592625, -0.34657381717123553];
                   composition = OptionsDict("T_e" => 0.5),
                   reactions=OptionsDict("charge_exchange_frequency"=>2*π*2.0 * sqrt(2)))

    # n_i=n_n T_e=4
    @long run_test(test_input_chebyshev_split_3_moments, 2*π*1.9919 * sqrt(2), -2*π*0.2491 * sqrt(2),
                   [-2.7725887202750434, -2.7725769429261002, -2.772570052243528,
                    -2.772568040068626, -2.7725678857317253, -2.7725680400686223];
                   composition = OptionsDict("T_e" => 4.0))
    # CX=2*π*2.0 case with T_e=4 is too hard to converge, so skip
end

function runtests()
    # Create a temporary directory for test output
    test_output_directory = get_MPI_tempdir()

    @testset "sound wave" verbose=use_verbose begin
        println("sound wave tests")

        @testset "finite difference" begin
            test_input_finite_difference["output"]["base_directory"] = test_output_directory
            run_test_set_finite_difference()

            test_input_finite_difference_split_1_moment["output"]["base_directory"] = test_output_directory
            @long run_test_set_finite_difference_split_1_moment()

            test_input_finite_difference_split_2_moments["output"]["base_directory"] = test_output_directory
            @long run_test_set_finite_difference_split_2_moments()

            test_input_finite_difference_split_3_moments["output"]["base_directory"] = test_output_directory
            run_test_set_finite_difference_split_3_moments()
        end

        @testset "Chebyshev" begin
            test_input_chebyshev["output"]["base_directory"] = test_output_directory
            run_test_set_chebyshev()

            test_input_chebyshev_split_1_moment["output"]["base_directory"] = test_output_directory
            run_test_set_chebyshev_split_1_moment()

            test_input_chebyshev_split_2_moments["output"]["base_directory"] = test_output_directory
            run_test_set_chebyshev_split_2_moments()

            test_input_chebyshev_split_3_moments["output"]["base_directory"] = test_output_directory
            run_test_set_chebyshev_split_3_moments()
        end
    end

    if global_rank[] == 0
        # Delete output directory to avoid using too much disk space
        rm(realpath(test_output_directory); recursive=true)
    end
end

end # SoundWaveTests


using .SoundWaveTests

SoundWaveTests.runtests()
