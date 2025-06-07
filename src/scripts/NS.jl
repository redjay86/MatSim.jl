cDir = @__DIR__
cd(cDir)
push!(LOAD_PATH,joinpath(cDir,"../libraries"))
#cd("C:\\Users\\rexja\\OneDrive - BYU-Idaho\\Classes\\Nelson Research Group\\MatSim.jl\\src\\libraries")
push!(LOAD_PATH,pwd())

using ase
using LennardJones
using BenchmarkTools
using YAML
using nSampling
# Load data in to test
resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")
input = YAML.load_file(joinpath(cDir,"NS.yml"))
LJ_average = LennardJones.get_LJ_averages(resultsPath)
myNS = nSampling.initialize(input["params"],LJ_average);# 10k allocations.  Need to optimize this.
nSampling.tune_step_sizes!(myNS,LJ_average)
nSampling.run_NS(myNS,LJ_average)