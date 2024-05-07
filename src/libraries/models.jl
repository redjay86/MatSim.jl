# Initialize directly from the YAML file.
function initializeModel(path::String)
    input = YAML.load_file(path)
    modelDict = Dict("lj"=>MatSim.initializeLJ)#,"sw"=>MatSim.initializeSW)
    return modelDict[lowercase(input["model"]["type"])](path)

end
