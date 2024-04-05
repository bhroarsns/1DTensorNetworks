using ITensors

function setupDir(target::String)
    resultdir = "./results/$(target)"
    snapshotdir = "./snapshots/$(target)"
    texdir = "./tex/$(target)"
    plotdir = "./plots/$(target)"
    if !isdir(resultdir)
        mkpath(resultdir)
    end
    if !isdir(snapshotdir)
        mkpath(snapshotdir)
    end
    if !isdir(texdir)
        mkpath(texdir)
    end
    if !isdir(plotdir)
        mkpath(plotdir)
    end
    return resultdir, snapshotdir
end

function snapshot(filename::String, tensor::ITensor, si::Index{Int}, lb::Index{Int}, rb::Index{Int})
    open(filename, "w") do io
        println(io, "# $(tags(si)), $(tags(lb)), $(tags(rb)), real, imag, abs, angle")
        for is in eachval(si)
            for il in eachval(lb)
                for ir in eachval(rb)
                    entry = tensor[si => is, lb => il, rb => ir]
                    println(io, "$(is), $(il), $(ir), $(real(entry)), $(imag(entry)), $(abs(entry)), $(angle(entry))")
                end
            end
        end
    end
    return nothing
end