cDir = @__DIR__
cd(cDir)
push!(LOAD_PATH,joinpath(cDir,"../libraries"))

using ase
using LennardJones
using nSampling
using YAML
using BenchmarkTools
using TimerOutputs

# Initialize the model
resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")
to = TimerOutput()
# Build model with averages of the draws for fit parameters.
LJ_average = LennardJones.get_LJ_averages(resultsPath)
input = YAML.load_file(joinpath(cDir,"NS.yml"))
myNS = nSampling.initialize(input["params"],LJ_average);# Initialize the simulation...
sEnergies =  reverse(sortperm([i.energies[2] for i in myNS.walkers]))
E_max = myNS.walkers[sEnergies[myNS.n_cull]].energies[2]
@btime nSampling.tune_step_sizes!($myNS,$LJ_average)
@btime nSampling.walk_single_walker!($myNS.walkers[5],$LJ_average,$myNS.walker_params,$E_max)
@btime ase.eval_energy($myNS.walkers[5],$LJ_average,force_recalc = true)
nSampling.do_cell_volume_step(myNS.walkers[5],LJ_average,myNS.walker_params,E_max,to)
@btime ase.DirectToCartesian!($myNS.walkers[5])
@btime ase.CartesianToDirect!($myNS.walkers[5])
show(to)
@time nSampling.propose_volume_step(myNS.walkers[5],myNS.walker_params.volume_step_size)
@btime ase.precalc_LJ_distances!($myNS.walkers[5],$LJ_average.cutoff)
@btime ase.get_neighbor_distances($myNS.walkers[5],6.0)
using StaticArrays
σ = SVector{3,Float64}(LJ_average.σ)
ϵ = SVector{3,Float64}(LJ_average.ϵ)
@btime vcat(-1.0 * $ϵ .* $σ.^6, $ϵ .* $σ.^12)
println(input)
nSampling.do_cell_shear_step(myNS.walkers[5],LJ_average,myNS.walker_params,E_max)
#buffers = nSampling.initialize_buffers()
ase.DirectToCartesian!(myNS.walkers[4])
println(myNS.walkers[4].positions)
ase.set_atoms_random!(myNS.walkers[3],1.5)
ase.atom_inside_cutoff(myNS.walkers[3],myNS.walkers[3].positions[1],2.0)
ase.get_neighbor_distances(myNS.walkers[8],3.5)
ase.precalc_LJ_distances!(myNS.walkers[8],LJ_average.cutoff)
println(myNS.walkers[3].lj_vec)
energy = ase.eval_PE_LJ(myNS.walkers[5],LJ_average,force_recalc = true)
println([x.energies[2] for x in myNS.walkers])
display(myNS.walkers[14].energies[2])
using BenchmarkTools
using StaticArrays
using Distributions
using LinearAlgebra

nAtoms = 8
check = SVector{nAtoms,Int64}([1 for x = 1:nAtoms])
@btime lowercase($myNS.walkers[5].coordSys[1][1]) == 'd'
println(myNS.walkers[5].lVecs)
using LinearAlgebra
mul!(myNS.walkers[5].lVecs, myNS.walkers[5].lVecs,myNS.walkers[5].latpar)
println(myNS.walker_params.shear_step_size)
@btime nSampling.do_cell_shear_step($myNS.walkers[5],$LJ_average,$myNS.walker_params,$myNS.walkers[5].model_energy * 10)
println(myNS.walkers[5].coordSys)
@btime ase.DirectToCartesian!($myNS.walkers[5])
@btime ase.CartesianToDirect!($myNS.walkers[5])
@code_warntype nSampling.propose_shear_step(myNS.walkers[5],myNS.walker_params.volume_step_size)
@btime nSampling.tune_step_sizes!($myNS,$LJ_average)
nSampling.run_NS(myNS,LJ_average)

using StaticArrays
check = false
pre = MMatrix{3,3,Float64,9}([1.0 0.0 1.2
       5.0 2.0 0.5
       2.1 4.2 1.4])
@allocated begin a = MMatrix{3,3,Float64,9}(undef)
 a .= pre
end
@allocated 
# Disregard all below
v1 = @view myNS.walkers[5].lVecs[:,2]

myNS.walkers[5].lVecs
display(ase.eval_KE(myNS.walkers[sEnergies[end-3]]))
println(myNS.walkers[sEnergies[end-3]].coordSys)
myNS.walkers[4].lj_vec
ase.eval_forces(myNS.walkers[1],LJ_average)
display(ase.get_neighbor_distances(myNS.walkers[1],4.0))
println(myNS.walker_params.MD_time_step)
sEnergies =  reverse(sortperm([i.energies[2] for i in myNS.walkers]))
E_max = myNS.walkers[sEnergies[myNS.n_cull]].energies[2]
nSampling.tune_step_sizes!(myNS.walkers[sEnergies[end-2]],LJ_average)

ase.eval_KE(myNS.walkers[4])
sum(0.5 * [x' * x for x in myNS.walkers[4].velocities])
nSampling.run_NS(myNS,LJ_average)
check =ase.get_nn_list(myNS.walkers[2],5.6)
check = myNS.walkers[1].positions
myvec = ["ab", "dc", "ef"]
!("ab" in myvec)
println(myNS.walkers[1].lVecs)
println(ase.cell_volume(myNS.walkers[1]))
ase.DirectToCartesian!(myNS.walkers[1])
myNS.walker_params.MD_time_step = 0.000001
ase.do_MD!(myNS.walkers[3],LJ_average,myNS.walker_params,E_max)
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
using StaticArrays
rand(SVector{3,Float64})


function check(orig_atoms)

       b = ase.atoms(orig_atoms.title,
             orig_atoms.latpar,
             orig_atoms.lVecs,
             orig_atoms.nType,
             orig_atoms.nAtoms,
             orig_atoms.coordSys,
             orig_atoms.positions,
             orig_atoms.velocities,
             orig_atoms.masses,
             orig_atoms.atomTypes,
             orig_atoms.species,
             orig_atoms.energies,
             orig_atoms.order,
             orig_atoms.lj_vec)
       b.positions[1] = SVector(100,100,100)
       return b
end
using StaticArrays
display(myNS.walkers[3].positions[1][1])
new_atoms = check(myNS.walkers[3])






function tester(a)
       b = MMatrix(a)
       b[:,1] .= [1.0, 2.2, 3.2]
       return SMatrix(b)

end

a = @SMatrix [5.0 2.3 9.2
              1.2 9.1 9.9
              8.3 7.5 2.6]
@btime tester($a)



possible = (do_cell_volume_step, do_cell_shear_step, do_cell_stretch_step)