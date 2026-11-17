module MomentKineticsDebugTests

include("setup.jl")

function runtests()
    @testset "moment_kinetics tests" begin
        include("kinetic_electron_tests.jl")
        include("sound_wave_tests.jl")
        include("fokker_planck_collisions_tests.jl")
        include("wall_bc_tests.jl")
        include("restart_interpolation_tests.jl")
        include("recycling_fraction_tests.jl")
        include("coll_krook_tests.jl")
        include("gyroaverage_tests.jl")

        manufactured_solutions_ext = Base.get_extension(moment_kinetics, :manufactured_solns_ext)
        if manufactured_solutions_ext !== nothing
            include("mms_tests.jl")
        end
    end
end

end # MomentKineticsDebugTests

using .MomentKineticsDebugTests

MomentKineticsDebugTests.runtests()
