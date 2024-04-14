using ITensors

function mkpathINE(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

function setupDir(target::String)
    resultdir = "./results/$(target)"
    snapshotdir = "./snapshots/$(target)"
    rm(snapshotdir, recursive=true)
    mkpathINE(resultdir)
    mkpathINE(snapshotdir)
    mkpathINE("./tex/$(target)")
    mkpathINE("./plots/$(target)")
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