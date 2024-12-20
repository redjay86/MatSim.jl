
#DataSet(dSet) = fromCrystals(dSet,standardize)

function fromdSet(dSet,standardize,fitTo)
    if standardize
        println("Standardizing Data")
        if lowercase(fitTo) == "peratom"
            println("Using per atom energies")
            for i in dSet.crystals
                i.fitEnergy = (i.energyPerAtomFP - dSet.meanEnergy)/dSet.stdEnergy-dSet.offset
            end
        elseif lowercase(fitTo) == "total"
            println("Using total energies")
            for i in dSet.crystals
                i.fitEnergy = (i.energyPerAtomFP * i.nAtoms - dSet.meanEnergy)/dSet.stdEnergy-dSet.offset
            end
        elseif lowercase(fitTo) == "fenth"
            println("Using formation energies")
            for i in dSet.crystals
                i.fitEnergy = (i.formationEnergyFP - dSet.meanEnergy)/dSet.stdEnergy-dSet.offset
            end
        end
        return DataSet(dSet.title,dSet.crystals,dSet.stdEnergy,dSet.meanEnergy,dSet.offset,dSet.nData,true,lowercase(fitTo))
    else
        println("Not standardizing data")
        if lowercase(fitTo) == "peratom"
            println("Using per atom energies")
            for i in dSet.crystals
                i.fitEnergy = i.energyPerAtomFP
            end
        elseif lowercase(fitTo) == "total"
            println("Using total energies")
            for i in dSet.crystals
                i.fitEnergy = i.energyPerAtomFP * i.nAtoms
            end
        elseif lowercase(fitTo) == "fenth"
            println("Using formation energies")
            for i in dSet.crystals
                i.fitEnergy = i.formationEnergyFP
            end
        end
        return DataSet(dSet.title,dSet.crystals,dSet.stdEnergy,dSet.meanEnergy,dSet.offset,dSet.nData,false,lowercase(fitTo))

    end

end
function getTraining_Holdout_Sets(dset::DataSet,nStructures,fitTo,standardize)

    training = sample(1:length(dset.crystals),nStructures,replace = false)
    holdout = setdiff(1:length(dset.crystals),training)

    if lowercase(fitTo) == "peratom"
        meanTraining = mean([i.energyPerAtomFP for i in dset.crystals[training]])
        stdTraining = std([i.energyPerAtomFP for i in dset.crystals[training]])

        meanHoldout = mean([i.energyPerAtomFP for i in dset.crystals[holdout]])
        stdHoldout = std([i.energyPerAtomFP for i in dset.crystals[holdout]])
    elseif lowercase(fitTo) == "total"
        meanTraining = mean([i.energyPerAtomFP * i.nAtoms for i in dset.crystals[training]])
        stdTraining = std([i.energyPerAtomFP * i.nAtoms for i in dset.crystals[training]])

        meanHoldout = mean([i.energyPerAtomFP * i.nAtoms for i in dset.crystals[holdout]])
        stdHoldout = std([i.energyPerAtomFP * i.nAtoms for i in dset.crystals[holdout]])
    elseif lowercase(fitTo) == "fenth"
        meanTraining = mean([i.formationEnergyFP for i in dset.crystals[training]])
        stdTraining = std([i.formationEnergyFP for i in dset.crystals[training]])

        meanHoldout = mean([i.formationEnergyFP for i in dset.crystals[holdout]])
        stdHoldout = std([i.formationEnergyFP for i in dset.crystals[holdout]])
    else
        error(" Can't tell what kind of energies you're wanting to fit to")
    end        


    trainingSet = DataSet(dset.title,dset.crystals[training],stdTraining,meanTraining,dset.offset,length(dset.crystals[training]),false,"?")
    holdoutSet = DataSet(dset.title,dset.crystals[holdout],stdHoldout,meanHoldout,dset.offset,length(dset.crystals[holdout]),false,"?")

    return fromdSet(trainingSet,standardize,fitTo), holdoutSet
    
end

function findPureEnergies(filePath)
    # Get the order of the system from the first POSCAR file.
    order = 0
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 8
            order = length(split(line))
            break
        end
        if idx > 100
            error("We shouldn't have to read this far into the structures.in file to find the order")

        end
    end
    println("order")
    println(order)
    pures = zeros(order)
    cLines = Vector{String}()
    counter = 0
    energyType = "peratom"
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 1
            energyType = split(line)[1]
            continue
        end
        if idx < 3
            continue
        end
        if occursin("#--",line)
            thisCrystal = Crystal(cLines,["N-a", "N-a"])
            energy = parse(Float64,cLines[end])
#            if lowercase(energyType) == "peratom"
#                thisCrystal.energyPerAtomFP = parse(Float64,cLines[end])
#                MatSim.formationEnergy!(thisCrystal,pures)
#            elseif lowercase(energyType) == "total"
#                thisCrystal.energyPerAtomFP = parse(Float64,cLines[end])/thisCrystal.nAtoms
#                MatSim.formationEnergy!(thisCrystal,pures)
#            elseif lowercase(energyType) == "fenth"
#                thisCrystal.formationEnergyFP = parse(Float64,cLines[end])
#                MatSim.totalEnergyFromFormationEnergy!(thisCrystal,pures)
#            else 
#                error("Can't recognize the first line of structures.in")
#            end
            if thisCrystal.nAtoms == 1
                pures[argmax(thisCrystal.nType)] = energy
            end
            cLines = Vector{String}()
            counter = 0
        else
            push!(cLines,line)
            counter += 1
        end
    end
    return pures
end
function readStructuresIn(filePath,species::Vector{String};overwriteLatPar = false,offset = 0.0)
    # Get the energies of the pure substances first so we can calculate formation energies
    pures = findPureEnergies(filePath)
    foundPures = true
    if all(pures .== 0)
        foundPures = false
        println(" Couldn't find pure energies. Formation Energies will not be accurate")
    end

    # The species list should always be in reverse alphabetical order. (All VASP calculations are performed under that convention)
    sort!(species,rev = true)

    data = Vector{Crystal}()
    cLines = Vector{String}()
    title = join(species,"-")
    counter = 0
    energyType = "peratom"
    for (idx,line) in enumerate(eachline(filePath))
        # The first line should tell me what kind of energies are present.
        if idx == 1
            energyType = split(line)[1]
            continue
        end
        if idx < 3
            continue
        end
        if occursin("#--",line)
            thisCrystal = Crystal(cLines,species,overwriteLatPar = overwriteLatPar)
            if lowercase(energyType) == "peratom"
                thisCrystal.energyPerAtomFP = parse(Float64,cLines[end])
                if foundPures
                    concentrations = thisCrystal.nType /thisCrystal.nAtoms
                    thisCrystal.formationEnergyFP = MatSim.formationEnergy(thisCrystal.energyPerAtomFP,pures,concentrations)
                else
                    thisCrystal.formationEnergyFP = NaN
                end
            elseif lowercase(energyType) == "total"
                thisCrystal.energyPerAtomFP = parse(Float64,cLines[end])/thisCrystal.nAtoms
                if foundPures
                    concentrations = thisCrystal.nType /thisCrystal.nAtoms
                    thisCrystal.formationEnergyFP = MatSim.formationEnergy(thisCrystal.energyPerAtomFP,pures,concentrations)
                else
                    thisCrystal.formationEnergyFP = NaN
                end
            elseif lowercase(energyType) == "fenth"
                thisCrystal.formationEnergyFP = parse(Float64,cLines[end])
                if foundPures
                    MatSim.totalEnergyFromFormationEnergy!(thisCrystal,pures)
                else
                    thisCrystal.energyPerAtomFP = NaN
                end
            else 
                error("Can't recognize the first line of structures.in")
            end
            if !isnan(thisCrystal.energyPerAtomFP)
                push!(data,thisCrystal)
            end
            cLines = Vector{String}()
            counter = 0
        else
            push!(cLines,line)
            counter += 1
        end
#        if occursin("#--",line)
#            nAtoms = sum(parse(Int64,x) for x in split(pos[idx + 6]))
#            startpoint = idx + 1
#            theend = idx + 7 + nAtoms
#            thisCrystal = Crystal(pos[startpoint:theend],species,overwriteLatPar = overwriteLatPar, energyPerAtomFP = parse(Float64,pos[theend + 2 ]))
#            if !isnan(thisCrystal.energyPerAtomFP)
#                push!(data,thisCrystal)
#            end
#        end
    end
    meanEnergy = mean([i.energyPerAtomFP for i in data])
    stdEnergy = std([i.energyPerAtomFP for i in data])    
    return DataSet(title,data,stdEnergy,meanEnergy,offset,length(data),false,"?")
end

function standardizeData!(dataSet::DataSet,offset::Float64)
#    meanEnergy = mean([i.energyPerAtomFP for i in dataSet.crystals])
#    stdEnergy = std([i.energyPerAtomFP for i in dataSet.crystals])
    #offset = 3
    for i in dataSet.crystals
        i.fitEnergy = (i.energyPerAtomFP - dataSet.meanEnergy)/dataSet.stdEnergy-offset
    end
#    dataSet.offset = offset
    dataSet.standardized = true
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
        write(io,string(crystal.energyPerAtomFP) * "\n")
    end
    close(io)
    
end


