using ITensors

function setupDir(target::String)
    resultdir = "./results/$(target)"
    snapshotdir = "./snapshots/$(target)"
    if !isdir(resultdir)
        mkpath(resultdir)
    end
    if !isdir(snapshotdir)
        mkpath(snapshotdir)
    end
    return resultdir, snapshotdir
end

function snapshot(filename::String, tensor::ITensor, siteInd::Index{Int}, leftbond::Index{Int}, rightbond::Index{Int})
    open(filename, "w") do io
        println(io, "# $(inds(siteInd)), $(inds(leftbond)), $(inds(rightbond)), real, imag, abs, angle")
        for isiteInd in eachval(siteInd)
            for ileftbond in eachval(leftbond)
                for irightbond in eachval(rightbond)
                    entry = tensor[siteInd => isiteInd, leftbond => ileftbond, rightbond => irightbond]
                    println(io, "$(isiteInd), $(ileftbond), $(irightbond), $(real(entry)), $(imag(entry)), $(abs(entry)), $(angle(entry))")
                end
            end
        end
    end
    return nothing
end