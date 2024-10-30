
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
using BenchmarkTools
using Profile
using Printf
cDir = @__DIR__
cd(cDir)
# Initialize the model
metrop,trainingSet,holdoutSet,LJ,dset = MatSim.initializeLJ(joinpath(cDir,"modelInputs-2.yml"));

display(dset.crystals[5].r6)
# Get samples
@time MatSim.getSamples(metrop,trainingSet,LJ)
path = joinpath(cDir,"struct_enum.out")

LJ_average = MatSim.LJAverages(metrop,LJ.order,LJ.cutoff)
display(LJ_average.ϵ)
MatSim.gss(path,LJ_average,dset.meanEnergy,dset.stdEnergy,holdoutSet.offset,readend = 20000)
# Display the results


MatSim.predPlot(metrop,LJ,trainingSet,holdoutSet,dset.meanEnergy,dset.stdEnergy,holdoutSet.offset)

MatSim.σ_hists(metrop,LJ)
hists = MatSim.std_hist(metrop,LJ)
MatSim.tracePlots(metrop,LJ)
display(plots)
MatSim.hists2d(metrop,LJ,"σ-ϵ";ar = :none)
display(plots[5])

using Plots
supertype(Plots.Plot)

using DelimitedFiles

data = readdlm(joinpath(cDir,"draws.out"),Float64,skipstart = 2)
data |> size|>typeof


using Printf
cDir = @__DIR__
cd(cDir)
var = 5.5
io = open(joinpath(cDir,"test.out"),"w")
@time printString = @sprintf "%15.5f" var
@time printString = repeat(" ",8) * string(var) * "\n"
@time write(io, repeat(" ",8) * string(var) * "\n")
@time @printf(io,"%15.5f",var)
@time Printf.format(Printf.Format("%15.5f"), var)
@time a = "Does this allocate"