
#using Pkg
#Pkg.installed()
#Pkg.activate("/Users/legoses/OneDrive - BYU-Idaho/codes/MatSim")
#Pkg.activate()
#Pkg.add("Plots")
#Pkg.resolve()
#Pkg.status()
#Pkg.add("YAML")
using Revise
using MatSim
#using Statistics
#using Distributions
#using Plots
#using YAML
cDir = @__DIR__
#cDir = pwd()
# Initialize the model
thisMetrop,trainingSet,holdoutSet,energyModel = MatSim.initializeModel(joinpath(cDir,"modelInputs.yml"))

# Run Metropolis
@time results = MatSim.getSamples(thisMetrop,trainingSet,energyModel)

# Look at the acceptance rates.
display(results.μ_accept)
println(results.σ_accept)

# Display the results
offset = 3.0
plots = MatSim.displayResults(results,energyModel,trainingSet,holdoutSet,offset = offset)

display(plots[1])
