include("../modules/iMPS.jl")
include("../modules/util.jl")

function doTEBD(
    # Hamiltonian
    modelname::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    sitetype::String,
    # TEBD parameters
    Δτs::Vector{Tuple{Float64,Int}},
    # MPS parameters
    D::Int,
    seed::Int,
    mpslen::Int=length(originalinds)
)
    target = "$(modelname)/iTEBD/mpslen=$(mpslen)/D=$(D)/seed=$(seed)"
    resultdir, snapshotdir = setupDir(target)
    open("$(snapshotdir)/spectrum.dat", "w")

    mps = randomInfiniteMPS(target, sitetype, D; mpslen, seed)
    normalize!(mps)
    takeSnapshot(mps, 0)
    β = 0.0
    totsteps = 0

    open("$(resultdir)/energy.dat", "w") do io
        println(io, "# D=$(D), seed=$(seed)")
        evs = real.(expectedvalues(mps, hloc, originalinds; normalized=true))
        println(io, "0, $(β), ", sum(evs) / length(evs), ", ", join(evs, ", "))

        for (Δτ, steps) in Δτs
            gate = exp(-Δτ * hloc)
            for istep in 1:steps
                print("\r", istep)
                update!(mps, gate, originalinds)
                normalize!(mps)
                takeSnapshot(mps, istep)
                evs = real.(expectedvalues(mps, hloc, originalinds; normalized=true))
                println(io, "$(totsteps+istep), $(β+Δτ*istep), ", sum(evs) / length(evs), ", ", join(evs, ", "))
            end

            β += Δτ * steps
            totsteps += steps
        end
    end
end