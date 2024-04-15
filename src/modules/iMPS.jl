using ITensors
using LinearAlgebra
using Random

include("util.jl")

mutable struct InfiniteMPS
    length::Int
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

function InfiniteMPS(siteTensors::Vector{ITensor}, bondWeights::Vector{ITensor})
    if length(siteTensors) != length(bondWeights)
        error("Length of siteTensors and bondWeights must be the same.")
    end
    mps = InfiniteMPS(length(siteTensors), siteTensors, bondWeights)
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

function takeSnapshot(mps::InfiniteMPS; nopr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing)
    for isite in 1:mps.length
        let stIO = ssio(nopr, "st:$(sitename(mps, isite))")
            if !isnothing(stIO)
                st = mps.siteTensors[isite]
                si = siteInd(mps, isite)
                bl, br = bondInds(mps, isite)
                println(stIO, "# $(tags(si)), $(tags(bl)), $(tags(br)), real, imag, abs, angle")
                for isi in eachval(si)
                    for ibl in eachval(bl)
                        for ibr in eachval(br)
                            entry = st[si=>isi, bl=>ibl, br=>ibr]
                            println(stIO, isi, ", ", ibl, ", ", ibr, ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
                        end
                    end
                end
            end
        end
    end

    for ibond in 1:mps.length
        let bwIO = ssio(nopr, "bw:$(bondname(mps, ibond))")
            if !isnothing(bwIO)
                bl, _, bw = bond(mps, ibond)
                println(bwIO, "# $(bondname(mps, ibond)), real, imag, abs, angle")
                for ibl in eachval(bl)
                    entry = bw[ibl, ibl]
                    println(bwIO, ibl, ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
                end
            end
        end
    end
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

function transfermatrix(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing)
    nopr = merge(opr, (methodcall="$(opr.methodcall)transfermatrix,",))
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

function symmetric(fp::ITensor)
    return (fp + swapprime(fp, 0, 1)) / 2.0
end

function antisymmetric(fp::ITensor)
    return (fp - swapprime(fp, 0, 1)) / 2.0
end

function symprojector(ind1::Index{Int}, ind2::Index{Int})
    D = dim(ind1)
    if dim(ind2) != D
        error("Given indices must have the same dimension.")
    end
    indsym = addtags(Index(D * (D + 1) ÷ 2, tags(ind1)), "Sym")
    indasym = addtags(Index(D * (D - 1) ÷ 2, tags(ind1)), "AntiSym")
    P = ITensor(ind1, ind2, indsym)
    Pinv = ITensor(ind1, ind2, indsym')
    Q = ITensor(ind1, ind2, indasym)
    Qinv = ITensor(ind1, ind2, indasym')

    for i in 0:D*(D+1)÷2-1
        qu = i ÷ (D + 1)
        re = rem(i, D + 1)
        qu2 = re ÷ (D - qu)
        lef = mod(re + 1, 1:D-qu)
        rig = lef + qu2 * D + (-1)^qu2 * (qu + qu2)
        if lef == rig
            P[ind1=>lef, ind2=>lef, indsym=>i+1] = 1.0
            Pinv[ind1=>lef, ind2=>lef, indsym'=>i+1] = 1.0
        else
            P[ind1=>lef, ind2=>rig, indsym=>i+1] = 1.0
            P[ind1=>rig, ind2=>lef, indsym=>i+1] = 1.0
            Pinv[ind1=>lef, ind2=>rig, indsym'=>i+1] = 0.5
            Pinv[ind1=>rig, ind2=>lef, indsym'=>i+1] = 0.5
            Q[ind1=>lef, ind2=>rig, indasym=>i+1-D] = 1.0
            Q[ind1=>rig, ind2=>lef, indasym=>i+1-D] = -1.0
            Qinv[ind1=>lef, ind2=>rig, indasym'=>i+1-D] = 0.5
            Qinv[ind1=>rig, ind2=>lef, indasym'=>i+1-D] = -0.5
        end
    end

    return P, Pinv, Q, Qinv, indsym, indasym
end

function sortSpectrum(D::ITensor, P::ITensor, spec::Spectrum, ep::Index{Int}, e::Index{Int})
    eigs = spec.eigs
    indord = sortperm(eigs; by=eig -> abs(eig), rev=true)
    spectrum = map(ie -> D[ie, ie], indord)
    λ = spectrum[begin]
    FPs = map(ie -> P * onehot(e => indord[ie]), findall(isapprox(λ), spectrum))
    return FPs, spectrum, λ
end

function fixedpoint(tm::ITensor, linkind::Index{Int}; calcAsym=true, opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing)
    nopr = merge(opr, (methodcall="$(opr.methodcall)fixedpoint,",))
    inds1 = [linkind, linkind']
    inds2 = uniqueinds(tm, inds1)

    PS, PSinv, PA, PAinv, indsym, indasym = symprojector(linkind, linkind')
    symtm = PS * real(tm) * replaceinds(PSinv, inds1, inds2)
    sFPs, sspec, sλ = sortSpectrum(eigen(symtm)...)
    sFPs = map(vec -> vec * PS, sFPs)
    degenFP = sFPs

    let sspecIO = ssio(nopr, "sspec")
        if !isnothing(sspecIO)
            for γ in sspec
                print(sspecIO, real(γ), ", ", imag(γ), ", ")
            end
            print(sspecIO, "\n")
            flush(sspecIO)
        end
    end

    if calcAsym
        asymtm = PA * real(tm) * replaceinds(PAinv, inds1, inds2)
        aFPs, aspec, aλ = sortSpectrum(eigen(asymtm)...)
        aFPs = map(vec -> vec * PA, aFPs)

        let errIO = ssio(nopr, "errtm")
            if !isnothing(errIO)
                symtmorg = replaceind(replaceind(PSinv, indsym', indsym) * symtm, indsym', indsym) * replaceinds(PS, inds1, inds2)
                asymtmorg = replaceind(replaceind(PAinv, indasym', indasym) * asymtm, indasym', indasym) * replaceinds(PA, inds1, inds2)
                errtm = norm(tm - symtmorg - asymtmorg) / norm(tm)
                println(errIO, errtm)
                flush(errIO)
            end
        end

        let aspecIO = ssio(nopr, "aspec")
            if !isnothing(aspecIO)
                for γ in aspec
                    print(aspecIO, real(γ), ", ", imag(γ), ", ")
                end
                print(aspecIO, "\n")
                flush(aspecIO)
            end
        end

        if sλ ≈ aλ
            append!(degenFP, aFPs)
        end
    end

    v = real(degenFP[begin])
    return v * sign(tr(v)), sλ, degenFP
end

function environment(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing)
    nopr = merge(opr, (methodcall="$(opr.methodcall)environment,",))
    El, Er, ll, lr = transfermatrix(mps, bondnum1, bondnum2; opr=nopr, ssio)
    σ, λl, degenFPl = fixedpoint(El, ll; opr=merge(nopr, (side="left",)), ssio)
    μ, λr, degenFPr = fixedpoint(Er, lr; opr=merge(nopr, (side="right",)), ssio)

    return σ, μ, ll, lr, (λl + λr) / 2.0, degenFPl, degenFPr
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

function canonicalize!(mps::InfiniteMPS, bondnum::Int; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing, fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    nopr = merge(opr, (methodcall="$(opr.methodcall)canonicalize!,",))
    mpslen = mps.length
    σ, μ, ll, lr, λ, degenFPl, degenFPr = environment(mps, bondnum; opr=nopr, ssio)
    Dl, Ul, _, el, _ = eigen(σ; ishermitian=true, cutoff=fpcutoff)
    Dr, Ur, _, _, _ = eigen(μ; ishermitian=true, cutoff=fpcutoff)

    let errσIO = ssio(nopr, "errσ")
        if !isnothing(errσIO)
            errσ = norm(σ - dag(Ul) * Dl * prime(Ul))
            println(errσIO, errσ)
            flush(errσIO)
        end
    end

    let errμIO = ssio(nopr, "errμ")
        if !isnothing(errμIO)
            errμ = norm(μ - dag(Ur) * Dr * prime(Ur))
            println(errμIO, errμ)
            flush(errμIO)
        end
    end

    let degenFPlIO = ssio(nopr, "degenFPl")
        if !isnothing(degenFPlIO)
            for fp in degenFPl
                println(degenFPlIO, fp)
                println(degenFPlIO, prime(dag(Ul)) * fp * Ul)
            end
            flush(degenFPlIO)
        end
    end

    let degenFPrIO = ssio(nopr, "degenFPr")
        if !isnothing(degenFPrIO)
            for fp in degenFPr
                println(degenFPrIO, fp)
                println(degenFPrIO, prime(dag(Ur)) * fp * Ur)
            end
            flush(degenFPrIO)
        end
    end

    bl, br, bw = bond(mps, bondnum)
    bn = bondname(mps, bondnum)

    Θ = sqrt.(Dl) * dag(Ul) * replaceinds(bw, [bl, br], [ll, lr]) * dag(Ur) * sqrt.(Dr)
    U, Σ, V = svd(Θ, el; cutoff=svcutoff, lefttags="Bond,$(bn),$(bn[1])", righttags="Bond,$(bn),$(bn[2])")

    let errΘIO = ssio(nopr, "errΘ")
        if !isnothing(errΘIO)
            errΘ = norm(Θ - U * Σ * V)
            println(errΘIO, errΘ)
            flush(errΘIO)
        end
    end

    X = Ul * inv.(sqrt.(Dl)) * U
    Y = V * inv.(sqrt.(Dr)) * Ur
    left = siteTensor(mps, bondnum) * replaceind(X, ll, bl)
    right = replaceind(Y, lr, br) * siteTensor(mps, bondnum + 1)

    mps.bondWeights[mod(bondnum, 1:mpslen)] = Σ
    mps.siteTensors[mod(bondnum, 1:mpslen)] = left
    mps.siteTensors[mod(bondnum + 1, 1:mpslen)] = right

    takeSnapshot(mps; nopr, ssio)

    return λ
end

function canonicalizeAll!(mps::InfiniteMPS; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing, fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    nopr = merge(opr, (methodcall="$(opr.methodcall)canonicalizeAll!,",))
    λs = Vector{Complex}(undef, mps.length)
    for ibond in eachindex(λs)
        λ = canonicalize!(mps, ibond; opr=merge(nopr, (bond=bondname(mps, ibond),)), ssio, fpcutoff, svcutoff)
        λs[ibond] = λ
    end
    return λs
end

function normalize!(mps::InfiniteMPS; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing, fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    nopr = merge(opr, (methodcall="$(opr.methodcall)normalize!,",))
    λs = canonicalizeAll!(mps; opr=nopr, ssio, fpcutoff, svcutoff)
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

    takeSnapshot(mps; nopr, ssio)

    return nothing
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}, firstsite::Int; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing, svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    nopr = merge(opr, (methodcall="$(opr.methodcall)update!,",))
    mpslen = mps.length
    gatelen = length(originalinds)
    if gatelen > mpslen
        error("Gate length must be smaller than MPS length.")
    end

    firstsite = mod(firstsite, 1:mpslen)
    # single site gate evolution do not require lbwinv / rbwinv
    if gatelen == 1
        Θ = mps.siteTensors[firstsite] * replaceind(gate, prime(originalinds[begin]), siteInd(mps, firstsite))
        mps.siteTensors[firstsite] = replaceind(Θ, originalinds[begin], siteInd(mps, firstsite))
        return 0.0
    end

    lastsite = firstsite + gatelen - 1
    minket, lbw, rbw, lbl, _ = contractKet(mps, firstsite, lastsite; minketonly=true)
    ketinds = uniqueinds(minket, lbw, rbw)

    tmp = sim(lbl)
    lbwinv = replaceind(inv.(lbw), lbl, tmp)
    rbwinv = inv.(rbw)
    Θ = replaceind(lbw, lbl, tmp) * minket * replaceinds(gate, prime.(originalinds), ketinds) * rbw
    Θ = replaceinds(Θ, originalinds, ketinds)
    Θorg = Θ

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
    mps.siteTensors[mod(lastsite, 1:mpslen)] = Θ

    let errUIO = ssio(nopr, "errU")
        if !isnothing(errUIO)
            Θnew, _ = contractKet(mps, firstsite, lastsite; minketonly=true)
            errU = norm(Θorg - Θnew) / norm(Θorg)
            println(errUIO, errU)
        end
    end

    mps.siteTensors[mod(lastsite, 1:mpslen)] = mps.siteTensors[mod(lastsite, 1:mpslen)] * rbwinv
    mps.siteTensors[firstsite] = lbwinv * mps.siteTensors[firstsite]
    takeSnapshot(mps; nopr, ssio)

    return nothing
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}; opr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> nothing, svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    nopr = merge(opr, (methodcall="$(opr.methodcall)update!,",))
    for firstsite in 1:mps.length
        update!(mps, gate, originalinds, firstsite; opr=nopr, ssio, svcutoff, newbonddim)
    end

    takeSnapshot(mps; nopr, ssio)
    return nothing
end

function randomInfiniteMPS(sitetype::String, bonddim::Int, mpslen::Int=2; seed::Union{Int,Nothing}=nothing)
    return randomInfiniteMPS(sitetype, ones(Int, mpslen) .* bonddim; seed)
end

function randomInfiniteMPS(sitetype::String, bonddims::Vector{Int}; seed::Union{Int,Nothing}=nothing)
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
    return InfiniteMPS(siteTensors, bondWeights)
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
