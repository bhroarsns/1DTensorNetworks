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

function symProjector(ind1::Index{Int}, ind2::Index{Int})
    D = dim(ind1)
    if dim(ind2) != D
        error("Given indices must have the same dimension.")
    end
    indsym = addtags(Index(D * (D + 1) ÷ 2, tags(ind1)), "Sym")
    P = ITensor(ind1, ind2, indsym)
    Pinv = ITensor(ind1, ind2, indsym')

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
        end
    end

    return P, Pinv
end

function getSymSector(tm::ITensor, linkind::Index{Int})
    inds1 = [linkind, linkind']
    inds2 = uniqueinds(tm, inds1)

    PS, PSinv = symProjector(linkind, linkind')
    return PS * real(tm) * replaceinds(PSinv, inds1, inds2), PS
end

function separator(ind1::Index{Int}, ind2::Index{Int})
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

function separateSector(tm::ITensor, linkind::Index{Int}; opr::Dict{String,String}=Dict{String,String}())
    inds1 = [linkind, linkind']
    inds2 = uniqueinds(tm, inds1)

    PS, PSinv, PA, PAinv, indsym, indasym = separator(linkind, linkind')
    symtm = PS * real(tm) * replaceinds(PSinv, inds1, inds2)
    asymtm = PA * real(tm) * replaceinds(PAinv, inds1, inds2)
    let (errIO, prefix) = ssio(opr, "errtm")
        if !isnothing(errIO)
            symtmorg = replaceind(PSinv, indsym', indsym) * symtm * replaceinds(PS, [indsym, inds1...], [indsym', inds2...])
            asymtmorg = replaceind(PAinv, indasym', indasym) * asymtm * replaceinds(PA, [indasym, inds1...], [indasym', inds2...])
            println(errIO, prefix..., norm(tm - symtmorg - asymtmorg) / norm(tm))
            flush(errIO)
        end
    end

    return symtm, PS, asymtm, PA
end

function peripheralSpectrum(tm::ITensor; findRealEV=true, toponly=false, opr::Dict{String,String}=Dict{String,String}())
    D, P, _, _, e = eigen(tm)
    sp = sortperm(storage(D); by=abs, rev=true)
    let (specIO, prefix) = ssio(opr, "spec")
        if !isnothing(specIO)
            print(specIO, prefix...)
            foreach(γ -> print(specIO, real(γ), ", ", imag(γ), ", "), storage(D)[sp])
            println(specIO)
            flush(specIO)
        end
    end
    spectralRadius = abs(storage(D)[sp[begin]])
    numdegen = count(isapprox(spectralRadius), abs.(storage(D)))
    pinds = sp[1:numdegen]
    pvecs = map(ie -> onehot(e => ie) * P, pinds)

    if findRealEV
        λind = @something findfirst(isreal, storage(D)[pinds]) error("No real eigenvalue found in the peripheral spectrum.")
        v = pvecs[λind]

        if toponly
            return real(storage(D)[pinds[λind]]), v
        else
            return spectralRadius, pvecs, real(storage(D)[pinds[λind]]), v
        end
    else
        return spectralRadius, pvecs
    end
end

function fixedpoint(tm::ITensor, linkind::Index{Int}; decomposed=true, opr::Dict{String,String}=Dict{String,String}(), plevel=0b1)
    nopr = mergewith(*, opr, Dict("methodcall" => "fixedpoint,"))

    if !decomposed
        symtm, PS = getSymSector(tm, linkind)
        λ, v = peripheralSpectrum(tm; toponly=true, opr=nopr)
        v = v * PS
        v = v / norm(v) * sign(tr(v))
        v, λ, tr(matrix(symtm)) / λ
    end

    symtm, PS, asymtm, PA = separateSector(tm, linkind)
    sρ, aρ, svects, avects, λ, v = if isodd(plevel)
        asymtask = @spawnat myid() + 1 peripheralSpectrum(asymtm; opr=merge(nopr, Dict("sector" => "asym")), findRealEV=false)
        sρ, svects, λ, v = peripheralSpectrum(symtm; opr=merge(nopr, Dict("sector" => "sym")))
        aρ, avects = fetch(asymtask)
        sρ, aρ, svects, avects, λ, v
    else
        sρ, svects, λ, v = peripheralSpectrum(symtm; opr=merge(nopr, Dict("sector" => "sym")))
        aρ, avects = peripheralSpectrum(asymtm; opr=merge(nopr, Dict("sector" => "asym")), findRealEV=false)
        sρ, aρ, svects, avects, λ, v
    end
    !(sρ ≈ aρ) && sρ < aρ && error("Spectral radius of the anti-symmetric sector of the TM exceeds that of the symmetric sector.")

    v = v * PS
    v = v / norm(v) * sign(tr(v))
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
            for sv in svects
                sFP = sv * PS
                sFP /= norm(sFP)
                println(degenFPIO, prime(dag(U)) * sFP * U)
            end
            if sρ ≈ aρ
                for av in avects
                    aFP = av * PA
                    aFP /= norm(aFP)
                    println(degenFPIO, prime(dag(U)) * aFP * U)
                end
            end
            flush(degenFPIO)
        end
    end

    return D, U, e, λ, tr(matrix(symtm)) / λ
end

function environment(mps::InfiniteMPS, bondnum1::Int, bondnum2::Int=bondnum1; decomposed=true, opr::Dict{String,String}=Dict{String,String}(), plevel::UInt8=0b11)
    nopr = mergewith(*, opr, Dict("methodcall" => "environment,"))
    El, Er, ll, lr = transfermatrix(mps, bondnum1, bondnum2; opr=nopr)
    if decomposed
        if isodd(plevel)
            lefttask = @spawnat myid() + 2^(count_ones(plevel >> 1)) fixedpoint(El, ll; decomposed, opr=merge(nopr, Dict("side" => "left")), plevel=plevel >> 1)
            Dr, Ur, _, λr, curTrR = fixedpoint(Er, lr; decomposed, opr=merge(nopr, Dict("side" => "right")), plevel=plevel >> 1)
            Dl, Ul, el, λl, curTrL = fetch(lefttask)
            return Dl, Ul, ll, Dr, Ur, lr, el, (λl + λr) / 2.0, (curTrL + curTrR) / 2.0
        else
            Dl, Ul, el, λl, curTrL = fixedpoint(El, ll; decomposed, opr=merge(nopr, Dict("side" => "left")), plevel=plevel >> 1)
            Dr, Ur, _, λr, curTrR = fixedpoint(Er, lr; decomposed, opr=merge(nopr, Dict("side" => "right")), plevel=plevel >> 1)
            return Dl, Ul, ll, Dr, Ur, lr, el, (λl + λr) / 2.0, (curTrL + curTrR) / 2.0
        end
    else
        if isodd(plevel)
            lefttask = @spawnat myid() + 1 fixedpoint(El, ll; decomposed, opr=merge(nopr, Dict("side" => "left")))
            μ, λr, curTrR = fixedpoint(Er, lr; decomposed, opr=merge(nopr, Dict("side" => "right")))
            σ, λl, curTrL = fetch(lefttask)
            return σ, ll, μ, lr, (λl + λr) / 2.0, (curTrL + curTrR) / 2.0
        else
            σ, λl, curTrL = fixedpoint(El, ll; decomposed, opr=merge(nopr, Dict("side" => "left")))
            μ, λr, curTrR = fixedpoint(Er, lr; decomposed, opr=merge(nopr, Dict("side" => "right")))
            return σ, ll, μ, lr, (λl + λr) / 2.0, (curTrL + curTrR) / 2.0
        end
    end
end
