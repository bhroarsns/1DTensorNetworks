function densitymatrix(mps::InfiniteMPS, dmlen::Int, firstsite::Int; normalized=false, opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "densitymatrix,"))
    lastsite = firstsite + dmlen - 1

    if normalized
        ket, minket, lbw, rbw, lb, _ = contractKet(mps, firstsite, lastsite; opr=nopr)
        ketinds = uniqueinds(ket, lbw, rbw)
        tmp = sim(lb)
        lbw = replaceind(lbw, lb, tmp)
        ket = lbw * ket * rbw

        # assume fixed points to be Kronecker's deltas
        rawDM = ket * prime(dag(ket); tags="Site")
        minket = if dmlen % mps.length == 0
            lbw
        else
            lbw * minket * rbw
        end
        bubble = minket * dag(minket)

        return rawDM / bubble[1], ketinds
    else
        repeatnum = dmlen ÷ mps.length
        ket, minket, lbw, rbw, lb, rb = contractKet(mps, firstsite, lastsite; opr=nopr)
        ketinds = uniqueinds(ket, lbw, rbw)

        # σ, μ, ll, lr, λ = environment(mps, firstsite - 1, lastsite)
        σ, ll, μ, lr, _, λ = environment(mps, firstsite - 1, lastsite; decomposed=false)
        lbw = replaceind(lbw, lb, ll)
        rbw = replaceind(rbw, rb, lr)

        σ2 = σ * lbw * prime(dag(lbw))
        μ2 = μ * rbw * prime(dag(rbw))
        rawDM = σ2 * ket * μ2 * prime(dag(ket))

        bubble = if dmlen % mps.length == 0
            minket = replaceind(lbw, rb, lr)
            σ * minket * μ * prime(dag(minket); tags="TMLink")
        else
            σ2 * minket * μ2 * prime(dag(minket); tags="Bond")
        end

        return rawDM / bubble[1] / λ^repeatnum, ketinds
    end
end

function expectedvalue(mps::InfiniteMPS, op::ITensor, originalinds::Vector{Index{Int}}, firstsite::Int; normalized=false, opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "expectedvalue,"))
    oplen = length(originalinds)
    ρ, ketinds = densitymatrix(mps, oplen, firstsite; normalized, opr=nopr)
    ev = ρ * replaceinds(op, [originalinds..., prime.(originalinds)...], [ketinds..., prime.(ketinds)...])
    return ev[1]
end

function expectedvalues(mps::InfiniteMPS, op::ITensor, originalinds::Vector{Index{Int}}; normalized=false, opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "expectedvalues,"))
    evs = Vector{Complex}(undef, mps.length)
    for firstsite in eachindex(evs)
        evs[firstsite] = expectedvalue(mps, op, originalinds, firstsite; normalized, opr=nopr)
    end
    return evs
end