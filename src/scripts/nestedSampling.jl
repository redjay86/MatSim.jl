using YAML
using MatSim
cDir = @__DIR__
cDir = pwd()
if !isfile(joinpath(cDir,"NS.yml"))
    error("Can't find the NS.yml file")
end
input = YAML.load_file(joinpath(cDir,"NS.yml"))
species = ["Ag", "Pt"]
myNS = MatSim.initializeSimulation(input["params"],species)# Initialize the simulation...

nestedSampling(myNS)