function transfermatrix(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1; opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "transfermatrix,"))
    bn1 = bondname(mps, bondnum1)
    bn2 = bondname(mps, bondnum2)
    minket1, lbw, rbw, lb, rb = contractKet(mps, bondnum1 + 1, bondnum1; minketonly=true, opr=nopr)
    minket2 = minket1
    if (bondnum2 - bondnum1) % mps.length != 0
        minket2, _, rbw, _, rb = contractKet(mps, bondnum2 + 1, bondnum2; minketonly=true, opr=nopr)
    end

    llink = replacetags(lb, "Bond", "TMLink")
    rlink = replacetags(rb, "Bond", "TMLink")

    lket = replaceind(lbw, lb, llink) * minket1
    rket = minket2 * replaceind(rbw, rb, rlink)
    El = lket * prime(dag(lket); tags=bn1)
    Er = rket * prime(dag(rket); tags=bn2)

    return El, Er, llink, rlink
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
            P[lef, lef, i+1] = 1.0
            Pinv[lef, lef, i+1] = 1.0
        else
            P[lef, rig, i+1] = 1.0
            P[rig, lef, i+1] = 1.0
            Pinv[lef, rig, i+1] = 0.5
            Pinv[rig, lef, i+1] = 0.5
            Q[lef, rig, i+1-D] = 1.0
            Q[rig, lef, i+1-D] = -1.0
            Qinv[lef, rig, i+1-D] = 0.5
            Qinv[rig, lef, i+1-D] = -0.5
        end
    end

    return P, Pinv, Q, Qinv, indsym, indasym
end

function sortSpectrum(tm::ITensor, proj::ITensor; notop=false, toponly=false, opr::Dict{String,String}=Dict{String,String}())
    D, P, spec, _, e = eigen(tm)
    eigs = spec.eigs
    indord = sortperm(eigs; by=eig -> abs(eig), rev=true)
    spectrum = storage(D)[indord]
    λind = indord[@something findfirst(isreal, spectrum) 1]
    λ = storage(D)[λind]

    if toponly
        v = onehot(e => λind) * P * proj
        return real(v / norm(v) * sign(tr(v))), λ
    end

    let (specIO, prefix) = ssio(opr, "spec")
        if !isnothing(specIO)
            print(specIO, prefix...)
            foreach(γ -> print(specIO, real(γ), ", ", imag(γ), ", "), spectrum)
            println(specIO)
            flush(specIO)
        end
    end

    FPs = map(ie -> begin
            FP = onehot(e => indord[ie]) * P * proj
            FP / norm(FP)
        end, findall(isapprox(abs(λ)), abs.(eigs)))

    if notop
        return FPs, λ
    else
        v = onehot(e => λind) * P * proj
        return FPs, real(v / norm(v) * sign(tr(v))), λ
    end
end

function fixedpoint(tm::ITensor, linkind::Index{Int}; decomposed=true, opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "fixedpoint,"))
    inds1 = [linkind, linkind']
    inds2 = uniqueinds(tm, inds1)

    PS, PSinv, PA, PAinv, indsym, indasym = symprojector(linkind, linkind')
    symtm = PS * real(tm) * replaceinds(PSinv, inds1, inds2)

    if !decomposed
        return sortSpectrum(symtm, PS; toponly=true, nopr)
    end

    asymtm = PA * real(tm) * replaceinds(PAinv, inds1, inds2)
    sFPs, v, sλ = sortSpectrum(symtm, PS; opr=merge(nopr, Dict("sector" => "sym")))
    aFPs, aλ = sortSpectrum(asymtm, PA; notop=true, opr=merge(nopr, Dict("sector" => "asym")))

    degenFP = if sλ ≈ aλ
        Iterators.flatten(zip(sFPs, aFPs))
    else
        sFPs
    end

    let (errIO, prefix) = ssio(nopr, "errtm")
        if !isnothing(errIO)
            symtmorg = replaceind(replaceind(PSinv, indsym', indsym) * symtm, indsym', indsym) * replaceinds(PS, inds1, inds2)
            asymtmorg = replaceind(replaceind(PAinv, indasym', indasym) * asymtm, indasym', indasym) * replaceinds(PA, inds1, inds2)
            println(errIO, prefix..., norm(tm - symtmorg - asymtmorg) / norm(tm))
            flush(errIO)
        end
    end

    D, U, _, e, _ = eigen(v; ishermitian=true)
    let (errvIO, prefix) = ssio(nopr, "errv")
        if !isnothing(errvIO)
            println(errvIO, prefix..., norm(v - dag(U) * D * prime(U)))
            flush(errvIO)
        end
    end
    let (degenFPIO, prefix) = ssio(nopr, "degenFP")
        if !isnothing(degenFPIO)
            print(degenFPIO, prefix...)
            foreach(fp -> begin
                    println(degenFPIO, fp)
                    println(degenFPIO, prime(dag(U)) * fp * U)
                end, degenFP)
            flush(degenFPIO)
        end
    end

    return D, U, e, sλ
end

function environment(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1; decomposed=true, opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "environment,"))
    El, Er, ll, lr = transfermatrix(mps, bondnum1, bondnum2; opr=nopr)
    if decomposed
        Dl, Ul, el, λl = fixedpoint(El, ll; decomposed, opr=merge(nopr, Dict("side" => "left")))
        Dr, Ur, _, λr = fixedpoint(Er, lr; decomposed, opr=merge(nopr, Dict("side" => "right")))
        return Dl, Ul, ll, Dr, Ur, lr, el, (λl + λr) / 2.0
    else
        σ, λl = fixedpoint(El, ll; decomposed, opr=merge(nopr, Dict("side" => "left")))
        μ, λr = fixedpoint(Er, lr; decomposed, opr=merge(nopr, Dict("side" => "right")))
        return σ, ll, μ, lr, (λl + λr) / 2.0
    end
end