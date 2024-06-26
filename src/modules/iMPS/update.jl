function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}, firstsite::Int; opr::Dict{String,String}=Dict{String,String}(), svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    nopr = mergewith(*, opr, Dict("methodcall" => "update!,"))
    mpslen = mps.length
    gatelen = length(originalinds)
    if gatelen > mpslen
        error("Gate length must be smaller than MPS length.")
    end

    firstsite = mod(firstsite, 1:mpslen)
    # single site gate evolution do not require lbwinv / rbwinv
    if gatelen == 1
        Θ = mps.siteTensors[firstsite] * replaceind(gate, originalinds[begin], siteInd(mps, firstsite))
        mps.siteTensors[firstsite] = replaceind(Θ, prime(originalinds[begin]), siteInd(mps, firstsite))
        return nothing
    end

    lastsite = firstsite + gatelen - 1
    minket, lbw, rbw, lbl, _ = contractKet(mps, firstsite, lastsite; minketonly=true, opr=nopr)
    ketinds = uniqueinds(minket, lbw, rbw)

    tmp = sim(lbl)
    lbwinv = replaceind(inv.(lbw), lbl, tmp)
    rbwinv = inv.(rbw)
    Θ = replaceind(lbw, lbl, tmp) * minket * replaceinds(gate, originalinds, ketinds) * rbw
    Θ = replaceinds(Θ, prime.(originalinds), ketinds)
    Θorg = Θ

    for ibond in firstsite:lastsite-1
        bn = bondname(mps, ibond)
        _, bi = bondInds(mps, ibond)
        U, _, V, spec, linku, linkv = svd(Θ, [tmp, ketinds[begin+ibond-firstsite]])

        let (specUIO, prefix) = ssio(merge(nopr, Dict("bond" => bn)), "uspec:$(bn)")
            if !isnothing(specUIO)
                print(specUIO, prefix...)
                foreach(γ -> print(specUIO, sqrt(γ), ", "), spec.eigs)
                println(specUIO)
                flush(specUIO)
            end
        end

        bonddim = isnothing(newbonddim) ? dim(bi) : min(newbonddim, dim(bi))
        truebonddim = bonddim
        for isv in bonddim+1:length(spec.eigs)
            if !(spec.eigs[isv] ≈ spec.eigs[bonddim])
                truebonddim = isv - 1
                break
            end
        end
        # print(", ", length(spec.eigs), "->", truebonddim, "(", bonddim, ")")
        bonddim = truebonddim
        bonddim = min(count(>(svcutoff^2), spec.eigs), bonddim)

        let (errUIO, prefix) = ssio(nopr, "errU")
            if !isnothing(errUIO)
                println(errUIO, prefix..., sum(spec.eigs[1:bonddim]), ", ", sum(spec.eigs))
                flush(errUIO)
            end
        end

        Σ = sqrt.(spec.eigs[1:bonddim])
        newu = Index(bonddim, "Bond,$(bn),$(bn[1])")
        tmp = Index(bonddim, "Bond,$(bn),$(bn[2])")
        mps.siteTensors[mod(ibond, 1:mpslen)] = U * δ(linku, newu)
        mps.bondWeights[mod(ibond, 1:mpslen)] = diagITensor(Σ, newu, tmp)
        Θ = δ(tmp, linkv) * V
    end
    mps.siteTensors[mod(lastsite, 1:mpslen)] = Θ * rbwinv
    mps.siteTensors[firstsite] = lbwinv * mps.siteTensors[firstsite]

    return nothing
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}; opr::Dict{String,String}=Dict{String,String}(), svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    nopr = mergewith(*, opr, Dict("methodcall" => "update!,"))
    for firstsite in 1:mps.length
        update!(mps, gate, originalinds, firstsite; opr=merge(nopr, Dict("fs" => "$(sitename(mps,firstsite))")), svcutoff, newbonddim)
    end

    return nothing
end