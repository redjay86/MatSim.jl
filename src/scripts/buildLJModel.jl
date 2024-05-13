
#using Pkg
#Pkg.installed()
#Pkg.activate("/Users/lancenelson/OneDrive - BYU-Idaho/codes/MatSim")
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
metrop,trainingSet,holdoutSet,LJ,dset = MatSim.initializeLJ(joinpath(cDir,"modelInputs-2.yml"))

# Run Metropolis

@time results = MatSim.getSamples(metrop,trainingSet,LJ)

# Look at the acceptance rates.
display(results.params_accept)
println(results.σ_accept)

# Display the results
plots = MatSim.displayResults(results,LJ,trainingSet,holdoutSet,dset.meanEnergy,dset.stdEnergy,holdoutSet.offset)

#gss()
display(plots[5])

