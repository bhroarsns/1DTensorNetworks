include("../modules/iMPS.jl")
include("../modules/util.jl")
using Printf

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

function recordErrs(mps::InfiniteMPS, errs::Vector{Vector{Float64}}, snapshotdir::String, errUs::Union{Vector{Float64},Nothing}=nothing)
    for ibond in eachindex(errs)
        open("./$(snapshotdir)/error/$(bondname(mps, ibond)).dat", "a") do io
            if isnothing(errUs)
                println(io, join(errs[ibond], ", "))
            else
                println(io, join(errs[ibond], ", "), ", ", errUs[ibond])
            end
        end
    end
end

function measurement(resultdir::String, mps::InfiniteMPS, hloc::ITensor, originalinds::Vector{Index{Int}}, istep::Int, β::Float64; singlesite::Union{ITensor,Nothing}=nothing, obs::Union{Vector{Tuple{ITensor,Vector{Index{Int}}}},Nothing}=nothing)
    evs = real.(expectedvalues(mps, hloc, originalinds; normalized=true))
    aveev = sum(evs) / length(evs)
    output = @sprintf("%03d, %+.16e, %s", istep, aveev, join(map(ev -> @sprintf("%+.16e", ev), evs), ", "))
    if !isnothing(singlesite)
        ssevs = real.(expectedvalues(mps, singlesite, [originalinds[begin]]; normalized=true))
        aveev += sum(ssevs) / length(ssevs)
        output = @sprintf("%03d, %+.16e, %s, %s", istep, aveev, join(map(ev -> @sprintf("%+.16e", ev), evs), ", "), join(map(ev -> @sprintf("%+.16e", ev), ssevs), ", "))
    end
    if !isnothing(obs)
        for (ob, orginds) in obs
            obevs = real.(expectedvalues(mps, ob, orginds; normalized=true))
            output *= @sprintf(", %+.16e, %s", sum(obevs) / length(obevs), join(map(ev -> @sprintf("%+.16e", ev), obevs), ", "))
        end
    end
    open("$(resultdir)/energy.dat", "a") do io
        println(io, output)
    end
end

function doiTEBD(
    # Hamiltonian
    modelname::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    sitetype::String,
    # TEBD parameters
    Δτs::Vector{Tuple{Float64,Int}},
    # MPS parameters
    D::Int,
    seed::Int;
    mpslen::Int=length(originalinds),
    singlesite::Union{ITensor,Nothing}=nothing,
    obs::Union{Vector{Tuple{ITensor,Vector{Index{Int}}}},Nothing}=nothing
)
    target = "$(modelname)/iTEBD/mpslen=$(mpslen)/D=$(D)/seed=$(seed)"
    resultdir, snapshotdir = setupDir(target)
    if !isdir("./$(snapshotdir)/spectrum")
        mkpath("./$(snapshotdir)/spectrum")
    end
    if !isdir("./$(snapshotdir)/error")
        mkpath("./$(snapshotdir)/error")
    end

    mps = randomInfiniteMPS(sitetype, D, mpslen; seed)
    open("$(resultdir)/energy.dat", "w") do io
        println(io, "# D=$(D), seed=$(seed)")
    end
    for ibond in 1:mpslen
        open("./$(snapshotdir)/spectrum/$(bondname(mps, ibond))_left.dat", "w")
        open("./$(snapshotdir)/spectrum/$(bondname(mps, ibond))_right.dat", "w")
        open("./$(snapshotdir)/error/$(bondname(mps, ibond)).dat", "w")
    end
    lspecs, rspecs, errs = normalize!(mps)
    recordSpecs(mps, lspecs, rspecs, snapshotdir)
    recordErrs(mps, errs, snapshotdir)
    takeSnapshot(mps, snapshotdir, 0)
    β = 0.0
    totsteps = 0
    measurement(resultdir, mps, hloc, originalinds, totsteps, β; singlesite, obs)

    for (Δτ, steps) in Δτs
        gate = exp(-Δτ * hloc)
        sgate = isnothing(singlesite) ? nothing : exp(-Δτ * singlesite)
        for istep in 1:steps
            print("\r", istep)
            errUs = update!(mps, gate, originalinds)
            lspecs, rspecs, errs = normalize!(mps)
            if !isnothing(sgate)
                update!(mps, sgate, [originalinds[begin]])
            end
            lspecs, rspecs, errs = normalize!(mps)
            recordSpecs(mps, lspecs, rspecs, snapshotdir)
            recordErrs(mps, errs, snapshotdir, errUs)
            takeSnapshot(mps, snapshotdir, istep)
            measurement(resultdir, mps, hloc, originalinds, totsteps + istep, β + Δτ * istep; singlesite, obs)
        end

        β += Δτ * steps
        totsteps += steps
    end

    return nothing
end