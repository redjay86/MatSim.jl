
function getTraining_Holdout_Sets(dset::DataSet,nStructures)

    training = sample(1:length(dset.crystals),nStructures,replace = false)
    holdout = setdiff(1:length(dset.crystals),training)

    trainingSet = DataSet(dset.crystals[training])
    holdoutSet = DataSet(dset.crystals[holdout])
    trainingSet, holdoutSet
end

function readStructuresIn(folder::String,file::String,species::Vector{String};overwriteLatPar = false)
    cd(folder)
    file = open(file,"r")
    pos = readlines(file)

    data = Vector{Crystal}()
    for (idx,line) in enumerate(pos)

        if occursin("#--",line)
            nAtoms = sum([parse(Int64,x) for x in split(pos[idx + 6])])
            startpoint = idx + 1
            theend = idx + 7 + nAtoms
            thisCrystal = Crystal(pos[startpoint:theend],species,overwriteLatPar = overwriteLatPar, energyFP = parse(Float64,pos[theend + 2 ]))
            if !isnan(thisCrystal.energyFP)
                push!(data,thisCrystal)
            end
        end
    end
    return DataSet(data)
end


function writeStructuresIn(path::String, structures::DataSet)
    io = open(path, "w")
    for crystal in structures.crystals
        write(io,"#-----------------------\n")
        write(io,crystal.title * "\n")
        write(io,string(crystal.latpar) * "\n")
        for lv in eachcol(crystal.lVecs)
            write(io,join(lv," ") * "\n")
        end
        write(io,join(crystal.nType," ") * "\n")
        write(io,crystal.coordSys[1] * "\n")
        for typ in crystal.atomicBasis
            for atom in typ
                write(io,join(atom," ") * "\n")
            end
        end
        write(io,"#Energy:\n")
        write(io,string(crystal.energyFP) * "\n")
    end
    close(io)
    
end