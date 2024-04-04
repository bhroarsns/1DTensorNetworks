using ITensors
using LinearAlgebra
using Random

include("util.jl")

mutable struct InfiniteMPS
    length::Int
    target::String
    siteTensors::Vector{ITensor}
    bondWeights::Vector{ITensor}
end

function siteTensor(mps::InfiniteMPS, sitenum::Int)
    return mps.siteTensors[mod(sitenum, 1:mps.length)]
end

function bondWeight(mps::InfiniteMPS, bondnum::Int)
    return mps.bondWeights[mod(bondnum, 1:mps.length)]
end

function sitename(mps::InfiniteMPS, sitenum::Int)
    return 'A' + mod(sitenum, 1:mps.length) - 1
end

function bondname(mps::InfiniteMPS, bondnum::Int)
    return string(sitename(mps, bondnum), sitename(mps, bondnum + 1))
end

function correctTags!(mps::InfiniteMPS)
    for site in eachindex(mps.siteTensors)
        st = siteTensor(mps, site)
        bwl = bondWeight(mps, site - 1)
        bwr = bondWeight(mps, site)

        bl = commonind(st, bwl)
        si = uniqueind(st, bwl, bwr)
        br = commonind(st, bwr)
        oldinds = [bl, si, br]

        blnew = removetags(bl, "Site,TMLink")
        blnew = addtags(blnew, string("Bond,", sitename(mps, site), ", ", bondname(mps, site - 1)))
        sinew = removetags(si, "Bond,TMLink")
        sinew = addtags(sinew, string("Site,", sitename(mps, site)))
        brnew = removetags(br, "Site,TMLink")
        brnew = addtags(brnew, string("Bond,", sitename(mps, site), ", ", bondname(mps, site)))
        newinds = [blnew, sinew, brnew]

        replaceinds!(st, oldinds, newinds)
        replaceinds!(bwl, oldinds, newinds)
        replaceinds!(bwr, oldinds, newinds)
    end
end

function InfiniteMPS(target::String, siteTensors::Vector{ITensor}, bondWeights::Vector{ITensor})
    if length(siteTensors) != length(bondWeights)
        error("Length of siteTensors and bondWeights must be the same.")
    end
    mps = InfiniteMPS(length(siteTensors), target, siteTensors, bondWeights)
    correctTags!(mps)
    return mps
end

function siteInd(mps::InfiniteMPS, sitenum::Int)
    return first(filter(hastags("Site"), inds(siteTensor(mps, sitenum))))
end

function bondInds(mps::InfiniteMPS, sitenum::Int)
    st = siteTensor(mps, sitenum)
    bwl = bondWeight(mps, sitenum - 1)
    bwr = bondWeight(mps, sitenum)
    return commonind(st, bwl), commonind(st, bwr)
end

function bond(mps::InfiniteMPS, bondnum::Int)
    bw = bondWeight(mps, bondnum)
    stl = siteTensor(mps, bondnum)
    str = siteTensor(mps, bondnum + 1)
    return commonind(bw, stl), commonind(bw, str), bw
end

function contractKet(mps::InfiniteMPS, firstsite::Int, lastsite::Int; minketonly=false)
    mpslen = mps.length
    ifirst = mod(firstsite, 1:mpslen)
    ilast = mod(lastsite, ifirst:ifirst+mpslen-1)
    repeatnum = (lastsite - firstsite) ÷ mpslen
    lbl, lbr, lbw = bond(mps, ifirst - 1)
    rbl, rbr, rbw = bond(mps, ilast)

    minket = deepcopy(siteTensor(mps, ifirst))
    for isite in ifirst+1:ilast
        minket = minket * bondWeight(mps, isite - 1) * siteTensor(mps, isite)
    end

    if minketonly
        return minket, lbw, rbw, lbl, rbr
    end

    # lbr - minket - rbl
    ket = deepcopy(minket)

    if repeatnum > 0
        # lbr - ket - lbl
        for isite in ilast+1:ifirst+mpslen-1
            ket = ket * bondWeight(mps, isite - 1) * siteTensor(mps, isite)
        end
        # lbr - ket - lbl * lbl - lbw - lbr'' -> lbr - ket - lbr''
        ket = ket * replaceind(lbw, lbr, lbr'')

        unitcell = deepcopy(ket)
        for irep in 1:repeatnum-1
            # lbr - ket - lbr'(2irep) * lbr'(2irep) - uc'(2irep) - lbr'(2irep+2) -> lbr - ket - lbr'(2irep+2)
            ket = ket * prime(unitcell, 2 * irep)
        end

        # lbr - ket - lbr'(2repeatnum) * lbr'(2repeatnum) - minket'(2repeatnum) - rbl'(2repeatnum)
        ket = ket * prime(minket, 2 * repeatnum)
        # lbr - ket - rbl
        ket = replaceind(ket, prime(rbl, 2 * repeatnum), rbl)
    end

    return ket, minket, lbw, rbw, lbl, rbr
end

function transfermatrix(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1)
    bn1 = bondname(mps, bondnum1)
    bn2 = bondname(mps, bondnum2)
    minket1, lbw, rbw, lb, rb = contractKet(mps, bondnum1 + 1, bondnum1; minketonly=true)
    minket2 = minket1
    if (bondnum2 - bondnum1) % mps.length != 0
        minket2, _, rbw, _, rb = contractKet(mps, bondnum2 + 1, bondnum2; minketonly=true)
    end

    llink = replacetags(lb, "Bond", "TMLink")
    rlink = replacetags(rb, "Bond", "TMLink")

    lket = replaceind(lbw, lb, llink) * minket1
    rket = minket2 * replaceind(rbw, rb, rlink)
    El = lket * prime(dag(lket); tags=bn1)
    Er = rket * prime(dag(rket); tags=bn2)
    return El, Er, llink, rlink
end

function fixedpoint(mps::InfiniteMPS, tm::ITensor, linkind::Index{Int})
    inds1 = [linkind, linkind']
    inds2 = uniqueinds(tm, inds1)

    D, P, _, _, e = eigen(tm, inds2, inds1)
    open("./snapshots/" * mps.target * "/spectrum.dat", "a") do io
        for ieig in eachval(e)
            print(io, real(D[ieig, ieig]), imag(D[ieig, ieig]))
        end
        print(io, "\n")
    end
    maxind = findlast(eig -> eig == maximum(eigs), eigs)
    λ = D[maxind, maxind]
    v = P * onehot(e => maxind)
    v = (v + swapprime(dag(v), 0, 1)) / 2.0
    return v * sign(tr(v)), λ
end

function environment(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1)
    El, Er, ll, lr = transfermatrix(mps, bondnum1, bondnum2)
    σ, λl = fixedpoint(mps, El, ll)
    μ, λr = fixedpoint(mps, Er, lr)

    return σ, μ, ll, lr, (λl + λr) / 2.0
end

function densitymatrix(mps::InfiniteMPS, dmlen::Int, firstsite::Int; normalized=false)
    lastsite = firstsite + dmlen - 1

    if normalized
        ket, minket, lbw, rbw, lb, _ = contractKet(mps, firstsite, lastsite)
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
        ket, minket, lbw, rbw, lb, rb = contractKet(mps, firstsite, lastsite)
        ketinds = uniqueinds(ket, lbw, rbw)

        σ, μ, ll, lr, λ = environment(mps, firstsite - 1, lastsite)
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

function expectedvalue(mps::InfiniteMPS, op::ITensor, originalinds::Vector{Index{Int}}, firstsite::Int; normalized=false)
    oplen = length(originalinds)
    ρ, ketinds = densitymatrix(mps, oplen, firstsite; normalized)
    ev = ρ * replaceinds(op, [prime.(originalinds)..., originalinds...], [ketinds..., prime.(ketinds)...])
    return ev[1]
end

function expectedvalues(mps::InfiniteMPS, op::ITensor, originalinds::Vector{Index{Int}}; normalized=false)
    evs = Vector{Complex}(undef, mps.length)
    for firstsite in eachindex(evs)
        evs[firstsite] = expectedvalue(mps, op, originalinds, firstsite; normalized)
    end
    return evs
end

function canonicalize!(mps::InfiniteMPS, bondnum::Int; fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    mpslen = mps.length
    σ, μ, ll, lr, λ = environment(mps, bondnum)
    Dl, Ul, _, el, _ = eigen(σ; ishermitian=true, cutoff=fpcutoff)
    Dr, Ur, _, _, _ = eigen(μ; ishermitian=true, cutoff=fpcutoff)

    bl, br, bw = bond(mps, bondnum)
    bn = bondname(mps, bondnum)

    Θ = sqrt.(Dl) * dag(Ul) * replaceinds(bw, [bl, br], [ll, lr]) * dag(Ur) * sqrt.(Dr)
    U, Σ, V = svd(Θ, el; cutoff=svcutoff, lefttags="Bond,$(bn),$(bn[1])", righttags="Bond,$(bn),$(bn[2])")

    X = Ul * inv.(sqrt.(Dl)) * U
    Y = V * inv.(sqrt.(Dr)) * Ur
    left = siteTensor(mps, bondnum) * replaceind(X, ll, bl)
    right = replaceind(Y, lr, br) * siteTensor(mps, bondnum + 1)

    mps.bondWeights[mod(bondnum, 1:mpslen)] = Σ
    mps.siteTensors[mod(bondnum, 1:mpslen)] = left
    mps.siteTensors[mod(bondnum + 1, 1:mpslen)] = right

    return λ
end

function canonicalizeAll!(mps::InfiniteMPS; fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    λs = Vector{Complex}(undef, mps.length)
    for ibond in eachindex(λs)
        λs[ibond] = canonicalize!(mps, ibond; fpcutoff, svcutoff)
    end
    return λs
end

function normalize!(mps::InfiniteMPS; fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    λs = canonicalizeAll!(mps; fpcutoff, svcutoff)
    divider = sqrt(sum(abs.(λs)) / length(λs))
    for ibond in eachindex(mps.bondWeights)
        bwnorm = norm(mps.bondWeights[ibond])
        divider /= bwnorm
        mps.bondWeights[ibond] /= bwnorm
    end
    divider ^= (1 / mps.length)

    for isite in eachindex(mps.siteTensors)
        mps.siteTensors[isite] /= divider
    end
    return nothing
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}, firstsite::Int; svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    mpslen = mps.length
    gatelen = length(originalinds)
    if gatelen > mpslen
        error("Gate length must be smaller than MPS length.")
    end

    firstsite = mod(firstsite, 1:mpslen)
    lastsite = firstsite + gatelen - 1

    minket, lbw, rbw, lbl, _ = contractKet(mps, firstsite, lastsite; minketonly=true)
    ketinds = uniqueinds(minket, lbw, rbw)
    tmp = sim(lbl)
    lbwinv = replaceind(inv.(lbw), lbl, tmp)
    rbwinv = inv.(rbw)
    Θ = replaceind(lbw, lbl, tmp) * minket * replaceinds(gate, prime.(originalinds), ketinds) * rbw
    Θ = replaceinds(Θ, originalinds, ketinds)

    for ibond in firstsite:lastsite-1
        bn = bondname(mps, ibond)
        _, bi = bondInds(mps, ibond)
        bonddim = isnothing(newbonddim) ? dim(bi) : newbonddim
        U, Σ, Θ, _, _, tmp = svd(Θ, [tmp, ketinds[begin+ibond-firstsite]];
            cutoff=svcutoff, maxdim=bonddim,
            lefttags="Bond,$(bn),$(bn[1])", righttags="Bond,$(bn),$(bn[2])")
        mps.siteTensors[mod(ibond, 1:mpslen)] = U
        mps.bondWeights[mod(ibond, 1:mpslen)] = Σ
    end

    mps.siteTensors[mod(lastsite, 1:mpslen)] = Θ * rbwinv
    mps.siteTensors[firstsite] = lbwinv * mps.siteTensors[firstsite]
    return nothing
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}; svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    for firstsite in 1:mps.length
        update!(mps, gate, originalinds, firstsite; svcutoff, newbonddim)
    end
    return nothing
end

function randomInfiniteMPS(target::String, sitetype, bonddim::Int, mpslen::Int=2; seed::Union{Int,Nothing}=nothing)
    return randomInfiniteMPS(target, sitetype, ones(Int, mpslen) .* bonddim; seed)
end

function randomInfiniteMPS(target::String, sitetype, bonddims::Vector{Int}; seed::Union{Int,Nothing}=nothing)
    mpslen = length(bonddims)
    siteTensors = Vector{ITensor}(undef, mpslen)
    bondWeights = Vector{ITensor}(undef, mpslen)
    bondindsl = Vector{Index{Int}}(undef, mpslen)
    bondindsr = Vector{Index{Int}}(undef, mpslen)
    for ibond in eachindex(bondWeights)
        bd = bonddims[ibond]
        lb = Index(bd)
        rb = Index(bd)
        bondWeights[ibond] = delta(lb, rb) / sqrt(bd)
        bondindsl[ibond] = lb
        bondindsr[ibond] = rb
    end
    if !isnothing(seed)
        Random.seed!(seed)
    end
    for isite in eachindex(siteTensors)
        si = siteind(sitetype)
        siteTensors[isite] = randomITensor(bondindsr[isite == 1 ? mpslen : isite - 1], si, bondindsl[isite])
    end
    return InfiniteMPS(target, siteTensors, bondWeights)
end

function takeSnapshot(mps::InfiniteMPS, number::Int)
    for isite in 1:mps.length
        if !isdir("./snapshots/$(mps.target)/" * sitename(mps, isite))
            mkpath("./snapshots/$(mps.target)/" * sitename(mps, isite))
        end
        st = mps.siteTensors[isite]
        si = siteInd(mps, isite)
        bl, br = bondInds(mps, isite)
        snapshot("./snapshots/$(mps.target)/" * sitename(mps, isite) * "/$(number).dat", st, si, bl, br)
    end
    for ibond in 1:mps.length
        if !isdir("./snapshots/$(mps.target)/" * bondname(mps, ibond))
            mkpath("./snapshots/$(mps.target)/" * bondname(mps, ibond))
        end
        open("./snapshots/$(mps.target)/" * bondname(mps, ibond) * "/$(number).dat", "w") do io
            bl, _, bw = bond(mps, ibond)
            for ibl in eachval(bl)
                entry = bw[ibl, ibl]
                println(io, "$(ibl), $(entry)")
            end
        end
    end
end


function getsitenum(mps::InfiniteMPS, sc::Char)
    sn = sc - 'A' + 1
    if 1 ≤ sn ≤ mps.length
        return sn
    else
        error("Invalid Site Character: \'$(sc)\'")
    end
end

function getbondnum(mps::InfiniteMPS, bn::String)
    leftsite = getsitenum(mps, bn[1])
    rightsite = getsitenum(mps, bn[2])
    if mod(leftsite + 1, 1:mps.length) == rightsite
        return leftsite
    else
        error("Invalid Bond Name: \"$(bondname)\"")
    end
end

function siteTensor(mps::InfiniteMPS, sc::Char)
    return siteTensor(mps, getsitenum(mps, sc))
end

function bondWeight(mps::InfiniteMPS, bn::String)
    return bondWeight(mps::InfiniteMPS, getbondnum(mps, bn))
end

function siteInd(mps::InfiniteMPS, sc::Char)
    return siteInd(mps, getsitenum(mps, sc))
end

function bondInds(mps::InfiniteMPS, sc::Char)
    return bondInds(mps, getsitenum(mps, sc))
end

function bond(mps::InfiniteMPS, bn::String)
    return bond(mps, getbondnum(mps, bn))
end

function contractKet(mps::InfiniteMPS, firstsite::Char, lastsite::Char; minketonly=false)
    return contractKet(mps, getsitenum(mps, firstsite), getsitenum(mps, lastsite); minketonly)
end

function transfermatrix(mps::InfiniteMPS, bn1::String, bn2::String=bn1)
    return transfermatrix(mps, getbondnum(mps, bn1), getbondnum(mps, bn2))
end

function environment(mps::InfiniteMPS, bn1::String, bn2::String=bn1)
    return environment(mps, getbondnum(mps, bn1), getbondnum(mps, bn2))
end

function densitymatrix(mps::InfiniteMPS, dmlen::Int, firstsite::Char; normalized=false)
    return densitymatrix(mps, dmlen, getsitenum(mps, firstsite); normalized)
end

function expectedvalue(mps::InfiniteMPS, op::ITensor, originalinds::Vector{Index{Int}}, firstsite::Char; normalized=false)
    return expectedvalue(mps, op, originalinds, getsitenum(mps, firstsite); normalized)
end

function canonicalize!(mps::InfiniteMPS, bn::String; fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    return canonicalize!(mps, getbondnum(mps, bn); fpcutoff, svcutoff)
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}, firstsite::Char; svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    return update!(mps, gate, originalinds, getsitenum(mps, firstsite); svcutoff, newbonddim)
end
