test_type = "Collisional Krook"
using moment_kinetics.type_definitions: OptionsDict

test_input = OptionsDict(
 "output" => OptionsDict(
     "run_name" => "coll_krook",
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
     "ngrid" => 1,
     "nelement" => 1,
    ),
 "z" => OptionsDict(
     "ngrid" => 4,
     "discretization" => "chebyshev_pseudospectral",
     "nelement" => 4,
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
     "ion_physics" => "coll_krook_ions",
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
     "nstep" => 2,
     "steady_state_residual" => true,
     "dt" => 1.0e-10,
     "nwrite" => 2,
     "type" => "SSPRK4",
     "nwrite_dfns" => 2,
     "print_nT_live" => true
    )
)

test_input_list = [test_input]
