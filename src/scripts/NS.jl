
cd("C:\\Users\\rexja\\OneDrive - BYU-Idaho\\Classes\\Nelson Research Group\\MatSim.jl\\src\\libraries")
push!(LOAD_PATH,pwd())

using ase
using PlotsMH
using LennardJones
using nSampling

using YAML
cDir = @__DIR__
cd(cDir)
# Initialize the model
resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")

# Plot the results
# Build model with averages of the draws for fit parameters.
LJ_average = PlotsMH.LJAverages(resultsPath)
#LennardJones.totalEnergy(LJ_average,["Ag","Pt"])
input = YAML.load_file(joinpath(cDir,"NS.yml"))
species = ["Ag", "Pt"]
fcc = ase.fromPOSCAR(joinpath(cDir,"POSCAR.fcc"),["Ag"])
sc = ase.fromPOSCAR(joinpath(cDir,"POSCAR.sc"),["Ag"])
myNS = nSampling.initialize(input["params"],species,LJ_average)# Initialize the simulation...
#walk_params = Dict("n_steps"=>myNS.n_single_walker_steps, "volume_step_size"=>myNS.volume_step_size, "shear_step_size" => myNS.shear_step_size, "stretch_step_size" => myNS.stretch_step_size,"max_volume_per_atom"=>myNS.max_volume_per_atom, "min_aspect_ratio"=>myNS.min_aspect_ratio)
check = myNS.walkers[1].positions
ase.DirectToCartesian!(myNS.walkers[1])
ase.do_MD!(myNS.walkers[1],20,LJ_average,10)
2 .* check #.* [2 * 5.0 ./ x for x in check] #.* 2 .* check
[2 * 5.0 ./ x for x in check] #.* 2 .* check
print(myNS.walkers[5].positions)
begin_energy = ase.eval_energy(myNS.walkers[2],LJ_average)
using Plots
using StaticArrays
using LinearAlgebra
nSampling.walk_single_walker!(myNS.walkers[2],LJ_average,myNS.walker_params,ase.eval_energy(myNS.walkers[2],LJ_average))
display(ase.eval_energy(myNS.walkers[2],LJ_average))
ase.gradientForce(LJ_average,myNS.walkers[2],@SVector[1,4],@SVector[2,2,2])
println(myNS.walkers[2].fitEnergy)
display(myNS.walkers[2].lVecs)
display(myNS.walkers[2].atomicBasis)
ase.mapIntoCell(myNS.walkers[2])
ase.DirectToCartesian!(myNS.walkers[2])

using LinearAlgebra
nHat = [[norm(x) > 0.0 ? x / norm(x) : x for x in y] for y in myNS.walkers[2].atomicBasis ]

any([isnan(z) for x in myNS.walkers[2].atomicBasis for y in x for z in y])
any(isnan.(myNS.walkers[2].atomicBasis[1]))
floor(-0.5)
rand(1:5)
#nSampling.GMC(myNS.configs[55],500,LJ_average)
nSampling.run_NS(myNS,LJ_average)
rands = nSampling.get_random_displacements([2,3])

println(rands[2])