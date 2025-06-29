module DataSets

using ase
using StatsBase
using Printf
using YAML
using QHull
using enumeration

struct DataSet
    title:: String
    configs::Vector{ase.atoms}
end


function initialize(file;species = ["N-A", "N-A"],offset = 0.0)
    # If you give a structures.in file directly, read directly from that. 
    if occursin("structures",file)
        println("Reading from" * file)
        dset = fromStructuresIn(file,species,overwriteLatPar = false,offset = offset)
        # Add other test cases as needed
    # Otherwise, parse the YAML file and read from the appropriate file.
    elseif occursin("yml",file)
        input = YAML.load_file(file,dicttype = Dict{String,Any})

        
        # Get the dataset
        species = String[x for x in input["dataset"]["species"]]
        dataFile = input["dataset"]["file"]
        println("Reading from" * dataFile)
        dset = fromStructuresIn(dataFile,species,overwriteLatPar = false,offset = offset)
    
    end
    return dset
end


function gss(file,model,species;readend = 100)
    pureatoms = ase.fccPures(species)
    pures = []
    for pure in pureatoms
        push!(pures,ase.eval_energy(pure,model))
    end
    enum=enumeration.read_header(file)

    cDir = dirname(file)
    io = open(joinpath(cDir,"gss.out"),"w")
    for idx in 1:readend
        atoms = ase.from_enum(file,idx,String.(species))
        atoms.model_energy = ase.eval_energy(atoms,model)
        conc = ase.get_concentrations(atoms)
        formationEnergyModel = ase.formationEnergy(atoms.model_energy,pures,conc)
        printString = @sprintf "%5d  %8.4f %8.4f %8.4f %8.4f\n" idx conc[1] conc[2] atoms.model_energy formationEnergyModel
        write(io,printString)
    end
    close(io)
end



function getConvexHull(file::String)

    nLines = countlines(file)  # Count the number of lines in the file.
    data = zeros(nLines,2)
    foundPureA = false
    foundPureB = false
    for (idx,line) in enumerate(eachline(file))  # Loop over each line in the file.
        concs = parse.(Float64,split(line)[2:3])  # Get the concentrations
        if isapprox(concs,[0.0,1.0])
            foundPureA = true
        elseif isapprox(concs,[1.0,0.0])
            foundPureB = true
        end
        data[idx,:] = parse.(Float64,split(line)[2:3:5])  # Get the second (concentration) and fifth (formation energy) entries in the file.
    #    push!(data,parse.(Float64,split(line)[2:3:5]))  # Get the second (concentration) and fifth (formation energy) entries in the file.
    end

    if !foundPureA
        data = vcat(data,[0.0 0.0])  # Add the pure A point if it is not present.
    end
    
    if !foundPureA
        data = vcat(data,[1.0 0.0])  # Add the pure A point if it is not present.
    end
    hull = QHull.chull(data)  # Find the convex hull of the data points.
    positiveRemoved = zeros(length(hull.vertices),2)  #Initialize an empty vector to hold the points that are below the zero formation energy line.
    counter = 1
    for (idx,value) in enumerate(hull.vertices) # Loop over the hull and only keep the points that are below the zero formation energy line.
        println(value)
        println(hull.points[value])
        if hull.points[value,2] <= 0.0 
            positiveRemoved[counter,:] = hull.points[value,:]  # Store the points that are below the zero formation energy line.
            counter += 1
        end
    end
    perms = sortperm(positiveRemoved[:,1])  # Sort the points by concentration.
    positiveRemoved = positiveRemoved[perms,:]  # Reorder the points by concentration.
    return unique(positiveRemoved,dims =1)  # Since I started with a matrix of zeros and then removed points, there may be unnecessary zeros remaining in the array. Remove them.
end

function undoRescaling!(dataSet::DataSet,fitTo)
    for i in dataSet.configs
        if fitTo == "peratom"
            i.energies[1] = i.energies[1] * i.nAtoms
        end
    end
end

# not used anymore, but kept for reference.
function standardizeData!(dataSet::DataSet)
    #meanEnergy = mean([i.FP_total_energy for i in dataSet.configs])
    #stdEnergy = std([i.FP_total_energy for i in dataSet.configs])
    #offset = 3
    for i in dataSet.configs
        i.energies[1] = (i.energies[1] - dataSet.meanEnergy)/dataSet.stdEnergy-dataSet.offset
    end
    #dataSet.offset = offset
#    dataSet.standardized = true
#    return nothing #stdEnergy, meanEnergy
end


function getTraining_Holdout_Sets(file::String;nTraining=100,species = ["N-A", "N-A"])
    # If you give a structures.in file directly, read directly from that. 
    if occursin("structures",file)
        println("Reading from" * file)
        dset = fromStructuresIn(file,species,overwriteLatPar = false)
        # Add other test cases as needed
    # Otherwise, parse the YAML file and read from the appropriate file.
    elseif occursin("yml",file)
        input = YAML.load_file(file,dicttype = Dict{String,Any})

        
        # Get the dataset
        species = String[x for x in input["dataset"]["species"]]
        dataFile = input["dataset"]["file"]
        offset = Float64(input["dataset"]["offset"])
        standardize = Bool(input["dataset"]["standardize"])
        fitTo = String(input["dataset"]["fitTo"])
        nTraining = Int(input["dataset"]["nTraining"])

        println("Reading from" * dataFile)
        dset = fromStructuresIn(dataFile,species,overwriteLatPar = false)

    end
    trainingSet,holdoutSet = getTraining_Holdout_Sets(dset,nTraining)

    return trainingSet,holdoutSet
end


function precalc_LJ_distances!(dset::DataSet,cutoff)
    for i in dset.configs
        ase.precalc_LJ_distances!(i,cutoff)
    end

end

function nConfigs(dset::DataSet)
    return length(dset.configs)
end
function getTraining_Holdout_Sets(dset::DataSet,nStructures)

    training = StatsBase.sample(1:length(dset.configs),nStructures,replace = false)
    holdout = setdiff(1:length(dset.configs),training)


    trainingSet = DataSet(dset.title,dset.configs[training])
    holdoutSet = DataSet(dset.title,dset.configs[holdout])

    return trainingSet, holdoutSet
    
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
    pures = Vector{ase.atoms}(undef,order)  # Initialize a vector to hold the pure substances.
    cLines = Vector{String}()
    counter = 0
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 1
            global energyType = split(line)[1]
            continue
        end
        if idx < 3
            continue
        end
        if occursin("#--",line)
            thisCrystal = ase.fromPOSCAR(cLines,["N-a", "N-a"])
            if lowercase(energyType) == "peratom"
                thisCrystal.energies[1]= parse(Float64,cLines[end]) * thisCrystal.nAtoms # Save the total energy.
            elseif lowercase(energyType) == "total"
                thisCrystal.energies[1] = parse(Float64,cLines[end])
            else 
                error("Can't recognize the first line of structures.in")
            end
            if count(!iszero,thisCrystal.nType) == 1  # If there is only one nonzero entry in the nType vector, it's a pure configuration.
                pures[argmax(thisCrystal.nType)] = thisCrystal
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
function fromStructuresIn(filePath,species::Vector{String};overwriteLatPar = false,offset = 0.0)

    # The species list should always be in reverse alphabetical order. (All VASP calculations are performed under that convention)
    sort!(species,rev = true)

    data = Vector{ase.atoms}()
    cLines = Vector{String}()
    title = join(species,"-")
    counter = 0
    for (idx,line) in enumerate(eachline(filePath))
        # The first line should tell me what kind of energies are present.
        if idx == 1
            global energyType = split(line)[1]
            continue
        end
        if idx < 3
            continue
        end
        if occursin("#--",line)
            thisCrystal = ase.fromPOSCAR(cLines,species,overwriteLatPar = overwriteLatPar)
            if lowercase(energyType) == "peratom"
                thisCrystal.energies[1] = parse(Float64,cLines[end]) * thisCrystal.nAtoms  # Always store the total energy, not per atom energy. 
            elseif lowercase(energyType) == "total"
                thisCrystal.FP_total_energy = parse(Float64,cLines[end])
            else 
                error("Can't recognize the first line of structures.in")
            end
            if !isnan(thisCrystal.energies[1])
                push!(data,thisCrystal)
            end
            cLines = Vector{String}()
            counter = 0
        else
            push!(cLines,line)
            counter += 1
        end

    end
    return DataSet(title,data)
end

function meanEnergy(dset::DataSet)
    return mean([i.FP_total_energy for i in dset.configs])
end

function stdEnergy(dset::DataSet)
    return std([i.FP_total_energy for i in dset.configs])
end
function rescaleData!(dataSet::DataSet,fitTo::String)
    # This function now just sets the energy value to the prescribed value before the fitting process.
    # The data set no longer gets modified using the mean and standard deviation.  That all happens in the model now.
    for i in dataSet.configs
        if fitTo == "peratom"
            i.energies[1] = i.energies[1]/i.nAtoms
        end
    end
end

function writeStructuresIn(path::String, structures::DataSet)
    io = open(path, "w")
    for crystal in structures.configs
        write(io,"#-----------------------\n")
        write(io,crystal.title * "\n")
        write(io,string(crystal.latpar) * "\n")
        for lv in eachcol(crystal.lVecs)
            write(io,join(lv," ") * "\n")
        end
        write(io,join(ase.nType(crystal)," ") * "\n")
        write(io,crystal.coordSys[1] * "\n")
        for typ in crystal.atomicBasis
            for atom in typ
                write(io,join(atom," ") * "\n")
            end
        end
        write(io,"#Energy:\n")
        write(io,string(crystal.FP_total_energy) * "\n")
    end
    close(io)
    
end



end