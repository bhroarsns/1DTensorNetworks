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

function countDegeneracy(eigs::Vector{Float64}, maxdim::Int)
    degen = Vector{Int}([1])
    for isv in Iterators.drop(eachindex(eigs), 1)
        if eigs[isv-1] ≈ eigs[isv]
            degen[end] += 1
        else
            push!(degen, 1)
        end
    end
    accdegen = filter(dim -> dim ≤ maxdim , accumulate(+, degen))
    newdim = isempty(accdegen) ? maxdim : last(accdegen)
    degen = collect(Iterators.flatten(map(num -> ones(Int, num) .* num, degen)))
    return degen, newdim
end

function diagonalizeRDM(curstep::Int, vec::ITensor, leftinds::Vector{Index{Int}}, maxdim::Int)
    UL, _, UR, spec, tmpl, tmpr = svd(vec, leftinds)
    degen, newdim = countDegeneracy(spec.eigs, maxdim)
    nlb = Index(newdim, "Bond,Left,$(curstep)-$(curstep+1)")
    nrb = Index(newdim, "Bond,Right,$(curstep)-$(curstep+1)")
    # truncation
    UL = UL * δ(tmpl, nlb)
    UR = UR * δ(tmpr, nrb)
    return UL, UR, spec.eigs, degen, nlb, nrb
end

function initDMRG(sitetype::String, hloc::ITensor, originalinds::Vector{Index{Int}}, maxdim::Int64; singlesite::Union{ITensor,Nothing}=nothing)
    plb = addtags(siteind(sitetype, 1), "Left")
    nls = addtags(siteind(sitetype, 2), "Left")
    nrs = addtags(siteind(sitetype, 2), "Right")
    prb = addtags(siteind(sitetype, 1), "Right")

    hlocinds = [prime.(originalinds)..., originalinds...]
    hL = replaceinds(hloc, hlocinds, [plb, nls, plb', nls'])
    hC = replaceinds(hloc, hlocinds, [nls, nrs, nls', nrs'])
    hR = replaceinds(hloc, hlocinds, [nrs, prb, nrs', prb'])

    if !isnothing(singlesite)
        orgind1 = originalinds[begin]
        hL += replaceinds(singlesite, [prime(orgind1), orgind1], [plb, plb']) * δ(nls, nls')
        hL += δ(plb, plb') * replaceinds(singlesite, [prime(orgind1), orgind1], [nls, nls'])
        hR += replaceinds(singlesite, [prime(orgind1), orgind1], [nrs, nrs']) * δ(prb, prb')
        hR += δ(nrs, nrs') * replaceinds(singlesite, [prime(orgind1), orgind1], [prb, prb'])
    end

    vec, gsEnergy = getGroundState(hL, hC, hR, [plb, nls, nrs, prb])
    UL, UR, eigs, degen, nlb, nrb = diagonalizeRDM(2, vec, [plb, nls], maxdim)

    return UL, UR, hL, hR, plb, nls, nlb, nrb, nrs, prb, eigs, degen, gsEnergy
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
    plb::Index,
    pls::Index,
    prs::Index,
    prb::Index,
    maxdim::Int;
    singlesite::Union{ITensor,Nothing}=nothing
)
    hL_new, nls = getNewBlock(curstep, sitetype, hloc, originalinds, false, hL, UL, pls, plb; singlesite)
    hR_new, nrs = getNewBlock(curstep, sitetype, hloc, originalinds, true, hR, UR, prs, prb; singlesite)
    hC_new = replaceinds(hloc, [prime.(originalinds)..., originalinds...], [nls, nrs, nls', nrs'])

    vec, gsEnergy = getGroundState(hL_new, hC_new, hR_new, [plb, nls, nrs, prb])
    UL_new, UR_new, eigs, degen, nlb, nrb = diagonalizeRDM(curstep, vec, [plb, nls], maxdim)

    return UL_new, UR_new, hL_new, hR_new, nls, nlb, nrb, nrs, eigs, degen, gsEnergy
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
