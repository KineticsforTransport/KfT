# provide option of running from command line via 'julia moment_kinetics.jl'
using Pkg
Pkg.activate(".")

using moment_kinetics
using moment_kinetics.type_definitions: OptionsDict

# Create a temporary directory for test output
test_output_directory = tempname()
mkpath(test_output_directory)

input = OptionsDict("output" => OptionsDict("run_name" => "precompilation",
                                            "base_directory" => test_output_directory),
                    "r" => OptionsDict("ngrid" => 5,
                                       "nelement" => 3,
                                       "discretization" => "chebyshev_pseudospectral"),
                    "inner_r_bc_1" => OptionsDict("bc" => "Neumann"),
                    "outer_r_bc_1" => OptionsDict("bc" => "Neumann"),
                    "z" => OptionsDict("ngrid" => 5,
                                       "nelement" => 3,
                                       "bc" => "wall",
                                       "discretization" => "chebyshev_pseudospectral"),
                    "vperp" => OptionsDict("ngrid" => 5,
                                           "nelement" => 3,
                                           "bc" => "zero",
                                           "L" => 4.0,
                                           "discretization" => "gausslegendre_pseudospectral"),
                    "vpa" => OptionsDict("ngrid" => 7,
                                         "nelement" => 3,
                                         "bc" => "zero",
                                         "L" => 8.0,
                                         "discretization" => "gausslegendre_pseudospectral"),
                    "vzeta" => OptionsDict("ngrid" => 5,
                                           "nelement" => 3,
                                           "bc" => "zero",
                                           "L" => 4.0,
                                           "discretization" => "chebyshev_pseudospectral"),
                    "vr" => OptionsDict("ngrid" => 5,
                                        "nelement" => 3,
                                        "bc" => "zero",
                                        "L" => 4.0,
                                        "discretization" => "chebyshev_pseudospectral"),
                    "vz" => OptionsDict("ngrid" => 7,
                                        "nelement" => 3,
                                        "bc" => "zero",
                                        "L" => 8.0,
                                        "discretization" => "chebyshev_pseudospectral"),
                    "composition" => OptionsDict("n_neutral_species" => 0),
                    "fokker_planck_collisions" => OptionsDict("use_fokker_planck" => true, "self_collisions" => true, "slowing_down_test" => true),
                    "timestepping" => OptionsDict("nstep" => 1, "dt" => 1.0e-11))

run_moment_kinetics(input)
