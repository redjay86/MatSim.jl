
using Pkg
#Pkg.installed()
#Pkg.activate("/Users/lancenelson/OneDrive - BYU-Idaho/codes/MatSim")
#Pkg.activate()
#Pkg.add("Plots")
#Pkg.resolve()
#Pkg.status()
#Pkg.add("YAML")
using Revise
using MatSim
using YAML
cDir = @__DIR__
cd(cDir)
# Initialize the model
metrop,trainingSet,holdoutSet,model,dset = MatSim.initializeLJ(joinpath(cDir,"modelInputs-2.yml"));

# Get samples
MatSim.getSamples(metrop,trainingSet,model)

resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")

# Plot the results
MatSim.predPlot(resultsPath,trainingSet)
MatSim.std_hist(resultsPath)
hists = MatSim.std_hist(resultsPath)
MatSim.tracePlots(resultsPath)
MatSim.hists2d(resultsPath,"σ-ϵ";ar = :none)

# Build model with averages of the draws for fit parameters.
LJ_average = MatSim.LJAverages(resultsPath)
path = joinpath(cDir,"struct_enum.out")
species = split(dset.title,"-")
MatSim.gss(path,LJ_average,species,readend = 2500)  # Predict energies of everything in struct_enum.out
MatSim.totalEnergy(dset.crystals[8],LJ_average)

if !isfile(joinpath(cDir,"NS.yml"))
    error("Can't find the NS.yml file")
end
input = YAML.load_file(joinpath(cDir,"NS.yml"))
species = ["Ag", "Pt"]
myNS = MatSim.initializeSimulation(input["params"],species,LJ_average)# Initialize the simulation...
println(myNS.configs[8].nAtoms)
MatSim.GMC(myNS.configs[55],500,LJ_average)
MatSim.nestedSampling(myNS,LJ_average)
display(dset.crystals[55].atomicBasis)
display(dset.crystals[55].coordSys)
display(LJ.meanEnergy)
display(model)
model = MatSim.LJ(model.order,model.cutoff, model.σ, model.ϵ,trainingSet.stdEnergy,trainingSet.meanEnergy, trainingSet.offset)

using QHull
using Pkg
using StaticArrays
convert(SVector{3}, 0.5 * rand(3) .- 0.25)
Pkg.add("QHull")


vaspDir = joinpath(cDir,"training_set")
sDir = joinpath(cDir,"structures.in")

MatSim.readVaspFolders(vaspDir,sDir,poscar = "POSCAR",energy = "total")

using YAML
using StaticArrays
using Plots
using LinearAlgebra
epsilon = 0.01
testCrystal = dset.crystals[3]
testCrystal.r6 .= 0.0
testCrystal.r12 .= 0.0
MatSim.DirectToCartesian!(testCrystal)
print(testCrystal.atomicBasis[2])
testCrystal.atomicBasis[2][1] -= SVector{3,Float64}(50 * epsilon,0,0)
loopBounds = SVector{3,Int64}(convert.(Int64,cld.(LJ_average.cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(testCrystal.latpar * testCrystal.lVecs)) )))

energies = []
forcesOne = []
forcesTwo = []
println(testCrystal.nType)
for i = 1:200
    energy = MatSim.totalEnergy(testCrystal,LJ_average)
    forceOne = MatSim.singleAtomForce(LJ_average,testCrystal,SVector(2,1),loopBounds)
    forceTwo = MatSim.gradientForce(LJ_average,testCrystal,[2,1],loopBounds)
    MatSim.DirectToCartesian!(testCrystal)
    testCrystal.atomicBasis[2][1] += SVector{3,Float64}(epsilon,0,0)
    push!(energies,energy)
    push!(forcesOne,forceOne[1])
    push!(forcesTwo,forceTwo[1])
end
length(energies)
x = -50 *epsilon:epsilon:149*epsilon
length(collect(x))
plot(x,energies)
plot(x,forcesOne)
plot(x,forcesTwo)

display(energies)


c1,c2 = MatSim.fccPures(["Ag","Pt"])
MatSim.writePOSCAR(c1,joinpath(cDir,"POSCAR.test"))
MatSim.writePOSCAR(c2,joinpath(cDir,"POSCAR2.test"))