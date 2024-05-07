
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
metrop,trainingSet,holdoutSet,LJ,dset = MatSim.initializeLJ(joinpath(cDir,"modelInputs.yml"))

LJ.params .= ones(2,2,2)
display(dset.crystals[8].ljvals)

check = MatSim.totalEnergy(dset.crystals[5],LJ)
# Run Metropolis
@time results = MatSim.getSamples(metrop,trainingSet,LJ)

# Look at the acceptance rates.
display(results.params_accept)
println(results.σ_accept)

# Display the results
plots = MatSim.displayResults(results,LJ,trainingSet,holdoutSet)

display(plots[2])


a = Dict("aa"=> Dict("eps"=>5.43, "σ"=> 2.8),"ab"=> Dict("eps"=>2.43, "σ"=> 3.8))
typeof(a)