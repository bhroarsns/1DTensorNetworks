include("../modules/iMPS.jl")
include("../modules/util.jl")
using Printf

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

function genssio(snapshotdir::String)
    return (opr::NamedTuple, ssname::String) -> begin
        mkpathINE("$(snapshotdir)/Step/$(opr.step)/$(opr.state)")
        mkpathINE("$(snapshotdir)/Spec/$(opr.state)")
        mkpathINE("$(snapshotdir)/Err/$(opr.state)")
        ssname == "errtm" && return open("$(snapshotdir)/Err/$(opr.state)/errtm.dat", "a"), (opr.step, ", ")
        ssname == "sspec" && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.bond)_$(opr.side)_sym.dat", "a"), (opr.step, ", ")
        ssname == "aspec" && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.bond)_$(opr.side)_asym.dat", "a"), (opr.step, ", ")
        ssname == "totspec" && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.bond)_$(opr.side)_tot.dat", "a"), (opr.step, ", ")
        ssname == "errσ" && return open("$(snapshotdir)/Err/$(opr.state)/errσ.dat", "a"), (opr.step, ", ")
        ssname == "errμ" && return open("$(snapshotdir)/Err/$(opr.state)/errμ.dat", "a"), (opr.step, ", ")
        ssname == "degenFPl" && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(opr.bond)_left_degenFP.dat", "w"), ""
        ssname == "degenFPr" && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(opr.bond)_right_degenFP.dat", "w"), ""
        ssname == "errΘ" && return open("$(snapshotdir)/Err/$(opr.state)/errΘ.dat", "a"), (opr.step, ", ")
        ssname == "errC" && return open("$(snapshotdir)/Err/$(opr.state)/errC.dat", "a"), (opr.step, ", ")
        ssname == "errU" && return open("$(snapshotdir)/Err/$(opr.state)/errU.dat", "a"), (opr.step, ", ", opr.fs, ", ")
        startswith(ssname, "uspec") && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.fs)_$(opr.bond).dat", "a"), (opr.step, ", ")
        if (opr.methodcall == "normalize!,") || (opr.methodcall == "update!,")
            startswith(ssname, "st") && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(ssname[4:end])", "w"), ""
            startswith(ssname, "bw") && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(ssname[4:end])", "w"), ""
        end
        return nothing, nothing
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
    open("$(resultdir)/energy.dat", "w") do io
        println(io, "# D=$(D), seed=$(seed)")
    end

    β = 0.0
    totsteps = 0
    mps = randomInfiniteMPS(sitetype, D, mpslen; seed)
    normalize!(mps; opr=(step=0, methodcall="",state="FUN"), ssio=genssio(snapshotdir))
    measurement(resultdir, mps, hloc, originalinds, totsteps, β; singlesite, obs)

    for (Δτ, steps) in Δτs
        gate = exp(-Δτ * hloc)
        sgate = isnothing(singlesite) ? nothing : exp(-Δτ * singlesite)
        for istep in 1:steps
            curstep = totsteps + istep
            print("\r", curstep)
            update!(mps, gate, originalinds; opr=(step=curstep, methodcall="",state="BSU"), ssio=genssio(snapshotdir))
            if !isnothing(sgate)
                normalize!(mps; opr=(step=curstep, methodcall="",state="BSN"), ssio=genssio(snapshotdir))
                update!(mps, sgate, [originalinds[begin]]; opr=(step=curstep, methodcall="",state="FUU"), ssio=genssio(snapshotdir))
            end
            normalize!(mps; opr=(step=curstep, methodcall="",state="FUN"), ssio=genssio(snapshotdir))
            measurement(resultdir, mps, hloc, originalinds, curstep, β + Δτ * curstep; singlesite, obs)
        end

        β += Δτ * steps
        totsteps += steps
    end

    println("\r", modelname, ": finished")

    return nothing
end