
function getTraining_Holdout_Sets(dset::DataSet,nStructures)

    training = sample(1:length(dset.crystals),nStructures,replace = false)
    holdout = setdiff(1:length(dset.crystals),training)

    meanTraining = mean([i.energyFP for i in dset.crystals[training]])
    stdTraining = std([i.energyFP for i in dset.crystals[training]])

    meanHoldout = mean([i.energyFP for i in dset.crystals[holdout]])
    stdHoldout = std([i.energyFP for i in dset.crystals[holdout]])



    trainingSet = DataSet(dset.crystals[training],stdTraining,meanTraining,dset.offset,length(dset.crystals[training]))
    holdoutSet = DataSet(dset.crystals[holdout],stdHoldout,meanHoldout,dset.offset,length(dset.crystals[holdout]))

#    if standardize
#        standardizeData!(trainingSet,offset)
#        standardizeData!(holdoutSet,offset)
#    end
    return trainingSet, holdoutSet
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
    meanEnergy = mean([i.energyFP for i in data])
    stdEnergy = std([i.energyFP for i in data])    
    return DataSet(data,stdEnergy,meanEnergy,0,length(data))
end

function standardizeData!(dataSet::DataSet,offset::Float64)
#    meanEnergy = mean([i.energyFP for i in dataSet.crystals])
#    stdEnergy = std([i.energyFP for i in dataSet.crystals])
    #offset = 3
    for i in dataSet.crystals
        i.energyFP = (i.energyFP - dataSet.meanEnergy)/dataSet.stdEnergy-offset
    end
    return nothing #stdEnergy, meanEnergy
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