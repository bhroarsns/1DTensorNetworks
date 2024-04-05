using ITensors
using LinearAlgebra

function initDMRG(hloc::ITensor, originalinds::Vector{Index{Int}}, sitetype::String, maxdim::Int64)
    plb = addtags(siteind(sitetype), "Left,n=1")
    nls = addtags(siteind(sitetype), "Left,n=2")
    nrs = addtags(siteind(sitetype), "Right,n=2")
    prb = addtags(siteind(sitetype), "Right,n=1")

    hL = replaceinds(hloc, [prime.(originalinds)..., originalinds...], [plb, nls, plb', nls'])
    hR = replaceinds(hloc, [prime.(originalinds)..., originalinds...], [nrs, prb, nrs', prb'])

    UL, UR, eigs, irreps, nlb, nrb, gsEnergy = getGroundState(
        2, hL,
        replaceinds(hloc, [prime.(originalinds)..., originalinds...], [nls, nrs, nls', nrs']),
        hR,
        plb, nls, nrs, prb, maxdim
    )

    return UL, UR, hL, hR, plb, nls, nlb, nrb, nrs, prb, eigs, irreps, gsEnergy
end

function getGroundState(
    curstep::Int,
    hL::ITensor,
    hC::ITensor,
    hR::ITensor,
    plb::Index,
    nls::Index,
    nrs::Index,
    prb::Index,
    maxdim::Int64
)
    hEff = hL * δ(nrs, nrs') * δ(prb, prb')
    hEff += δ(plb, plb') * hC * δ(prb, prb')
    hEff += δ(plb, plb') * δ(nls, nls') * hR

    _, U, spec, _, e = eigen(hEff; ishermitian=true)
    eigs = spec.eigs
    minind = findfirst(eig -> eig == minimum(eigs), eigs)
    gsEnergy = eigs[minind]
    vec = U * onehot(e => minind)
    UL, _, UR, spec, tl, tr = svd(vec, [plb, nls])
    eigs = spec.eigs
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
    nlb = Index(newdim, "Bond,Left,$(curstep)-$(curstep+1)")
    nrb = Index(newdim, "Bond,Right,$(curstep)-$(curstep+1)")
    UL = UL * δ(tl, nlb)
    UR = UR * δ(tr, nrb)
    irreps = collect(Iterators.flatten(map(num -> ones(Int, num) .* num, irreps)))
    return UL, UR, eigs, irreps, nlb, nrb, gsEnergy
end

function dmrgStep(
    curstep::Int,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    sitetype::String,
    UL::ITensor,
    UR::ITensor,
    hL::ITensor,
    hR::ITensor,
    plb::Index,
    pls::Index,
    prs::Index,
    prb::Index,
    maxdim::Int
)
    nls = addtags(siteind(sitetype), "Left,n=$(curstep)")
    nrs = addtags(siteind(sitetype), "Right,n=$(curstep)")

    hL_new = UL * hL * prime(dag(UL)) * delta(nls, nls')
    hL_new += UL * replaceinds(hloc, [prime.(originalinds)..., originalinds...], [pls, nls, pls', nls']) * replaceinds(dag(UL), [pls, plb], [pls', plb'])
    hR_new = UR * hR * prime(dag(UR)) * delta(nrs, nrs')
    hR_new += UR * replaceinds(hloc, [prime.(originalinds)..., originalinds...], [nrs, prs, nrs', prs']) * replaceinds(dag(UR), [prs, prb], [prs', prb'])
    hC_new = replaceinds(hloc, [prime.(originalinds)..., originalinds...], [nls, nrs, nls', nrs'])

    UL_new, UR_new, eigs, irreps, nlb, nrb, gsEnergy = getGroundState(curstep, hL_new, hC_new, hR_new, plb, nls, nrs, prb, maxdim)

    return UL_new, UR_new, hL_new, hR_new, nls, nlb, nrb, nrs, eigs, irreps, gsEnergy
end

