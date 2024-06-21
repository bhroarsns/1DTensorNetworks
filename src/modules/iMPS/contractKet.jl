function contractKet(mps::InfiniteMPS, firstsite::Int, lastsite::Int; minketonly=false, opr::Dict{String,String}=Dict{String,String}())
    mpslen = mps.length
    ifirst = mod(firstsite, 1:mpslen)
    ilast = mod(lastsite, ifirst:ifirst+mpslen-1)
    repeatnum = (lastsite - firstsite) ÷ mpslen
    lbl, lbr, lbw = bond(mps, ifirst - 1)
    rbl, rbr, rbw = bond(mps, ilast)

    minket :: ITensor = deepcopy(siteTensor(mps, ifirst))
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