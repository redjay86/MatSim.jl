

cd("/Users/legoses/OneDrive - BYU-Idaho/codes/MatSim/src/libraries/")
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
display(input)
species = ["Ag", "Pt"]
fcc = ase.fromPOSCAR(joinpath(cDir,"POSCAR.fcc"),["Ag"])
sc = atomsSim.fromPOSCAR(joinpath(cDir,"POSCAR.sc"),["Ag"])
myNS = nSampling.initialize(input["params"],species,LJ_average)# Initialize the simulation...
#walk_params = Dict("n_steps"=>myNS.n_single_walker_steps, "volume_step_size"=>myNS.volume_step_size, "shear_step_size" => myNS.shear_step_size, "stretch_step_size" => myNS.stretch_step_size,"max_volume_per_atom"=>myNS.max_volume_per_atom, "min_aspect_ratio"=>myNS.min_aspect_ratio)
begin_energy = ase.eval_energy(myNS.walkers[2],LJ_average)
nSampling.walk_single_walker!(myNS.walkers[2],LJ_average,myNS.walker_params,ase.eval_energy(myNS.walkers[2],LJ_average))
display(ase.eval_energy(myNS.walkers[2],LJ_average))
rand(1:5)
#nSampling.GMC(myNS.configs[55],500,LJ_average)
nSampling.run_NS(myNS,LJ_average)


