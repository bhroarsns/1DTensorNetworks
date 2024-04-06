include("../modules/iMPS.jl")
include("../modules/util.jl")

function takeSnapshot(mps::InfiniteMPS, snapshotdir::String, number::Int)
    for isite in 1:mps.length
        if !isdir("./$(snapshotdir)/" * sitename(mps, isite))
            mkpath("./$(snapshotdir)/" * sitename(mps, isite))
        end
        st = mps.siteTensors[isite]
        si = siteInd(mps, isite)
        bl, br = bondInds(mps, isite)
        snapshot("./$(snapshotdir)/" * sitename(mps, isite) * "/$(number).dat", st, si, bl, br)
    end
    for ibond in 1:mps.length
        if !isdir("./$(snapshotdir)/" * bondname(mps, ibond))
            mkpath("./$(snapshotdir)/" * bondname(mps, ibond))
        end
        open("./$(snapshotdir)/" * bondname(mps, ibond) * "/$(number).dat", "w") do io
            bl, _, bw = bond(mps, ibond)
            for ibl in eachval(bl)
                entry = bw[ibl, ibl]
                println(io, "$(ibl), $(entry)")
            end
        end
    end
end

function recordSpecs(mps::InfiniteMPS, lspecs::Vector{Vector{Complex}}, rspecs::Vector{Vector{Complex}}, snapshotdir::String)
    for ibond in eachindex(lspecs)
        open("./$(snapshotdir)/spectrum/$(bondname(mps, ibond))_left.dat", "a") do io
            println(io, join(map(λ -> "$(real(λ)), $(imag(λ))", lspecs[ibond]), ", "))
        end
        open("./$(snapshotdir)/spectrum/$(bondname(mps, ibond))_right.dat", "a") do io
            println(io, join(map(λ -> "$(real(λ)), $(imag(λ))", rspecs[ibond]), ", "))
        end
    end
end

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
    if !isdir("./$(snapshotdir)/spectrum")
        mkpath("./$(snapshotdir)/spectrum")
    end

    mps = randomInfiniteMPS(sitetype, D, mpslen; seed)
    for ibond in 1:mpslen
        open("./$(snapshotdir)/spectrum/$(bondname(mps, ibond))_left.dat", "w")
        open("./$(snapshotdir)/spectrum/$(bondname(mps, ibond))_right.dat", "w")
    end
    lspecs, rspecs = normalize!(mps)
    recordSpecs(mps, lspecs, rspecs, snapshotdir)
    takeSnapshot(mps, snapshotdir, 0)
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
                lspecs, rspecs = normalize!(mps)
                recordSpecs(mps, lspecs, rspecs, snapshotdir)
                takeSnapshot(mps, snapshotdir, istep)
                evs = real.(expectedvalues(mps, hloc, originalinds; normalized=true))
                println(io, "$(totsteps+istep), $(β+Δτ*istep), ", sum(evs) / length(evs), ", ", join(evs, ", "))
            end

            β += Δτ * steps
            totsteps += steps
        end
    end

    return nothing
end