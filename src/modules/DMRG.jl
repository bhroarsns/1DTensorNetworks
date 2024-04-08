using ITensors
using LinearAlgebra

function getGroundState(hL::ITensor, hC::ITensor, hR::ITensor, his::Vector{Index{Int}})
    hEff = hL * δ(his[3], his[3]') * δ(his[4], his[4]')
    hEff += δ(his[1], his[1]') * hC * δ(his[4], his[4]')
    hEff += δ(his[1], his[1]') * δ(his[2], his[2]') * hR

    _, U, spec, _, e = eigen(hEff; ishermitian=true)
    vec, gsEnergy = let eigs = spec.eigs
        minind = findfirst(eig -> eig == minimum(eigs), eigs)
        vec = U * onehot(e => minind)
        vec, eigs[minind]
    end
    return vec, gsEnergy
end

function countIrreps(eigs::Vector{Float64}, maxdim::Int)
    irreps = Vector{Int}([1])
    for isv in Iterators.drop(eachindex(eigs), 1)
        if eigs[isv-1] ≈ eigs[isv]
            irreps[end] += 1
        else
            push!(irreps, 1)
        end
    end
    newdim = reduce((i, j) -> i + j > maxdim ? i : i + j, irreps)
    if newdim == 0
        newdim = maxdim
    end
    irreps = collect(Iterators.flatten(map(num -> ones(Int, num) .* num, irreps)))
    return irreps, newdim
end

function diagonalizeRDM(curstep::Int, vec::ITensor, leftinds::Vector{Index{Int}}, maxdim::Int)
    UL, _, UR, spec, tmpl, tmpr = svd(vec, leftinds)
    irreps, newdim = countIrreps(spec.eigs, maxdim)
    nleftbond = Index(newdim, "Bond,Left,$(curstep)-$(curstep+1)")
    nrightbond = Index(newdim, "Bond,Right,$(curstep)-$(curstep+1)")
    # truncation
    UL = UL * δ(tmpl, nleftbond)
    UR = UR * δ(tmpr, nrightbond)
    return UL, UR, spec.eigs, irreps, nleftbond, nrightbond
end

function initDMRG(sitetype::String, hloc::ITensor, originalinds::Vector{Index{Int}}, maxdim::Int64; singlesite::Union{ITensor,Nothing}=nothing)
    pleftbond = addtags(siteind(sitetype, 1), "Left")
    nleftsite = addtags(siteind(sitetype, 2), "Left")
    nrightsite = addtags(siteind(sitetype, 2), "Right")
    prightbond = addtags(siteind(sitetype, 1), "Right")

    hlocinds = [prime.(originalinds)..., originalinds...]
    hL = replaceinds(hloc, hlocinds, [pleftbond, nleftsite, pleftbond', nleftsite'])
    hC = replaceinds(hloc, hlocinds, [nleftsite, nrightsite, nleftsite', nrightsite'])
    hR = replaceinds(hloc, hlocinds, [nrightsite, prightbond, nrightsite', prightbond'])

    if !isnothing(singlesite)
        orgind1 = originalinds[begin]
        hL += replaceinds(singlesite, [prime(orgind1), orgind1], [pleftbond, pleftbond']) * δ(nleftsite, nleftsite')
        hL += δ(pleftbond, pleftbond') * replaceinds(singlesite, [prime(orgind1), orgind1], [nleftsite, nleftsite'])
        hR += replaceinds(singlesite, [prime(orgind1), orgind1], [nrightsite, nrightsite']) * δ(prightbond, prightbond')
        hR += δ(nrightsite, nrightsite') * replaceinds(singlesite, [prime(orgind1), orgind1], [prightbond, prightbond'])
    end

    vec, gsEnergy = getGroundState(hL, hC, hR, [pleftbond, nleftsite, nrightsite, prightbond])
    UL, UR, eigs, irreps, nleftbond, nrightbond = diagonalizeRDM(2, vec, [pleftbond, nleftsite], maxdim)

    return UL, UR, hL, hR, pleftbond, nleftsite, nleftbond, nrightbond, nrightsite, prightbond, eigs, irreps, gsEnergy
end

function dmrgStep(
    curstep::Int,
    sitetype::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    UL::ITensor,
    UR::ITensor,
    hL::ITensor,
    hR::ITensor,
    pleftbond::Index,
    pleftsite::Index,
    prightsite::Index,
    prightbond::Index,
    maxdim::Int;
    singlesite::Union{ITensor,Nothing}=nothing
)
    hL_new, nleftsite = getNewBlock(curstep, sitetype, hloc, originalinds, false, hL, UL, pleftsite, pleftbond; singlesite)
    hR_new, nrightsite = getNewBlock(curstep, sitetype, hloc, originalinds, true, hR, UR, prightsite, prightbond; singlesite)
    hC_new = replaceinds(hloc, [prime.(originalinds)..., originalinds...], [nleftsite, nrightsite, nleftsite', nrightsite'])

    vec, gsEnergy = getGroundState(hL_new, hC_new, hR_new, [pleftbond, nleftsite, nrightsite, prightbond])
    UL_new, UR_new, eigs, irreps, nleftbond, nrightbond = diagonalizeRDM(curstep, vec, [pleftbond, nleftsite], maxdim)

    return UL_new, UR_new, hL_new, hR_new, nleftsite, nleftbond, nrightbond, nrightsite, eigs, irreps, gsEnergy
end

function getNewBlock(
    curstep::Int,
    sitetype::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    isright::Bool,
    hprev::ITensor,
    U::ITensor,
    ps::Index,
    pb::Index;
    singlesite::Union{ITensor,Nothing}=nothing
)
    ns = addtags(siteind(sitetype, curstep), (isright ? "Right" : "Left"))
    #  ┌───U───pb      pb───U───┐
    #  │   │                │   │
    #  #   ps              ps   #
    #  #   ps  ns      ns  ps   #
    #  │   │   │        │   │   │
    # ┌─────┐  │        │  ┌─────┐
    # │hprev│  │        │  │hprev│
    # └─────┘  │        │  └─────┘
    #  │   │   │        │   │   │
    #  #'  ps' ns'     ns' ps'  #'
    #  #'  ps'             ps'  #'
    #  │   │                │   │
    #  └───Ū───pb'     pb'──Ū───┘
    hnew = U * hprev * prime(dag(U)) * delta(ns, ns')

    sites = isright ? [ns, ps] : [ps, ns]
    nhl = replaceinds(hloc, [prime.(originalinds)..., originalinds...], [sites..., prime.(sites)...])
    #  ┌───U───pb      pb───U───┐
    #  │   │                │   │
    #  #   ps              ps   #
    #  │   ps  ns      ns  ps   │
    #  │   │   │        │   │   │
    #  │  ┌─────┐      ┌─────┐  │
    #  │  │ nhl │      │ nhl │  │
    #  │  └─────┘      └─────┘  │
    #  │   │   │        │   │   │
    #  │   ps' ns'     ns' ps'  │
    #  #   ps'             ps'  #
    #  │   │                │   │
    #  └───Ū───pb'     pb'──Ū───┘
    hnew += U * nhl * replaceinds(dag(U), [ps, pb], [ps', pb'])

    if !isnothing(singlesite)
        ss = replaceinds(singlesite, [prime(originalinds[begin]), originalinds[begin]], [ns, ns'])
        #  ┌───U───pb      pb───U───┐      
        #  │   │                │   │      
        #  #   ps              ps   #     
        #  │   │   ns      ns   │   │      pb   ns   ns   pb 
        #  │   │   │        │   │   │       │   │     │   │ 
        #  │   │ ┌──┐      ┌──┐ │   │       │ ┌──┐   ┌──┐ │ 
        #  │   │ │ss│      │ss│ │   │   =   │ │ss│   │ss│ │ 
        #  │   │ └──┘      └──┘ │   │       │ └──┘   └──┘ │ 
        #  │   │   │        │   │   │       │   │     │   │ 
        #  │   │   ns'     ns'  │   │      pb'  ns'  ns'  pb'
        #  #   ps              ps   #     
        #  │   │                │   │      
        #  └───Ū───pb'     pb'──Ū───┘      
        hnew += delta(pb, pb') * ss
    end

    return hnew, ns
end
