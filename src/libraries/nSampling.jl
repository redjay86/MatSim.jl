module nSampling

using ase
using LennardJones
using DataSets:DataSet
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using TimerOutputs


mutable struct NS_walker_params
    n_single_walker_steps:: Int64
    MC_atom_step_size:: Float64
    MD_time_step:: Float64
    MD_reject_eps:: Float64
    n_atom_steps:: Int64
    n_cell_volume_steps:: Int64
    n_cell_shear_steps:: Int64
    n_cell_stretch_steps:: Int64
    atom_traj_length:: Int64
    min_aspect_ratio:: Float64
    max_volume_per_atom:: Float64
    volume_step_size:: Float64
    shear_step_size:: Float64
    stretch_step_size:: Float64
    KE_max:: Float64
    atom_algorithm:: String
end

struct NS
    n_walkers:: Int
    n_cull:: Int
    n_iter:: Int
    cell_P:: Float64
    eps:: Float64
    walker_params:: NS_walker_params
    walkers::Vector{ase.atoms}
end

# I added this for performance enhancements. There are several places
# in the code where we calculate a new cell geometry and/or transformation matrices.  
# This happens a lot in NS so let's just initialize them once and modify them in place
# rather than creating them every single time.
function initialize_buffers()
    return NS_buffers(
        MMatrix{3,3,Float64,9}(zeros(3,3)),
        MMatrix{3,3,Float64,9}(zeros(3,3)),
        MMatrix{3,3,Float64,9}(zeros(3,3)))

end

function NS_walker_params(inputs)
    final_settings = Dict()
    for (key,value) in NS_defaults()
        if !(key in keys(inputs))
            # Use the default
            final_settings[key] = value
        else
            # Use custom setting
            final_settings[key] = inputs[key]
        end

    end
    nAtoms = sum(inputs["n_Atoms"])
    max_lc = (nAtoms * inputs["max_volume_per_atom"])^(1/3)  # This is the maximum lattice constant

     return NS_walker_params(final_settings["n_single_walker_steps"],max_lc * final_settings["MC_atom_step_size"],final_settings["MD_time_step"],final_settings["MD_reject_eps"],final_settings["n_atom_steps"],
        final_settings["n_cell_volume_steps"],final_settings["n_cell_shear_steps"],final_settings["n_cell_stretch_steps"],
        final_settings["atom_traj_length"],final_settings["min_aspect_ratio"],final_settings["max_volume_per_atom"],
        final_settings["volume_step_size"],final_settings["shear_step_size"],final_settings["stretch_step_size"],inputs["KE_max"],final_settings["atom_algorithm"])
end

function NS_defaults()
    defaults = Dict("min_atom_separation"=> 1.0,
                    "n_walkers" => 100, 
                    "n_cull" => 1, 
                    "atom_traj_length"=> 8, 
                    "eps"=> 1e-6,
                    "n_iter" => 2000,
                    "atom_algorithm"=> "MD", 
                    "MC_atom_step_size" => 0.5, 
                    "MD_time_step"=>0.1,
                    "MD_reject_eps"=>1e-2,
                    "n_single_walker_steps"=> 600, 
                    "n_atom_steps"=> 100, 
                    "n_cell_volume_steps"=> 5, 
                    "n_cell_shear_steps"=>8,
                    "n_cell_stretch_steps"=>8, 
                    "cell_P"=>0.0, 
                    "min_aspect_ratio" => 0.8, 
                    "max_volume_per_atom"=>32.0,
                    "volume_step_size" => 0.5, 
                    "shear_step_size"=>0.1, 
                    "stretch_step_size"=>0.1)

    return defaults
end



function initialize(inputs,model)
    nTypes = inputs["n_Atoms"]
    nAtoms = sum(inputs["n_Atoms"])
    max_lc = (nAtoms * inputs["max_volume_per_atom"])^(1/3)  # This is the maximum atom separation.
    latticeConstant = (nAtoms * inputs["max_volume_per_atom"] * rand()^(1/(nAtoms + 1)))^(1/3)
    lVecs = SMatrix{3,3,Float64,9}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
#    error("stopping")
    configs = [ase.initialize_cell_shape(latticeConstant,lVecs,nTypes) for i in 1:inputs["n_walkers"]]
    if haskey(inputs, "cell_P")
        inputs["KE_max"] = 3/2 * inputs["cell_P"] * inputs["max_volume_per_atom"] * nAtoms
    end
    walker_params = NS_walker_params(inputs)
#    walk_params = Dict("n_steps" => 10,
#                       "shear_step_size" => inputs["shear_step_size"],
#                       "stretch_step_size" => inputs["stretch_step_size"])
    # Now do a random walk on the cell shape and place atoms at random positions
    save_n_steps = walker_params.n_single_walker_steps
    walker_params.n_single_walker_steps = 10
    for iConfig = 1:inputs["n_walkers"]
        configs[iConfig] = do_cell_shape_walk!(configs[iConfig],walker_params)
        ase.set_atoms_random!(configs[iConfig],inputs["min_atom_separation"])
#        error("Stop")
        ase.set_masses(configs[iConfig],1.0)
        if lowercase(inputs["atom_algorithm"]) == "md"
            ase.set_random_unit_velocities!(configs[iConfig],walker_params.KE_max)
        else
            ase.set_velocities(configs[iConfig],0.0)
        end
        configs[iConfig].energies[2] = ase.eval_energy(configs[iConfig],model,force_recalc = true)
    end
    walker_params.n_single_walker_steps = save_n_steps
    return NS(inputs["n_walkers"],inputs["n_cull"],inputs["n_iter"],inputs["cell_P"],inputs["eps"],
              walker_params,configs)
end

function adjust_step_sizes!(walk_params,key,rate)

    if rate > 1//2
        if key == "volume"
            walk_params.volume_step_size *=1.05
        elseif key == "shear"
            walk_params.shear_step_size *=1.05
        elseif key == "stretch"
            walk_params.stretch_step_size *=1.05
        elseif key == "atoms"
            if walk_params.atom_algorithm == "MD"
                walk_params.MD_time_step *=1.02
            else
                walk_params.MC_atom_step_size *=1.05
            end
        end
        return false
    elseif rate < 1//4
        if key == "volume"
            walk_params.volume_step_size *=0.95
        elseif key == "shear"
            walk_params.shear_step_size *=0.95
        elseif key == "stretch"
            walk_params.stretch_step_size *=0.95
        elseif key == "atoms"
            if walk_params.atom_algorithm == "MD"
                walk_params.MD_time_step *=0.98
            else
                walk_params.MC_atom_step_size *=0.95
            end
        end
        return false
    else
        return true
    end

end

function tune_step_sizes!(NS,model::LennardJones.model)
    sEnergies =  reverse(sortperm([i.energies[2] for i in NS.walkers]))
    keeps = sEnergies[NS.n_cull + 1: end]
    E_max = NS.walkers[sEnergies[1]].energies[2]
    atoms = deepcopy(NS.walkers[sample(keeps,1)[1]])
   
    allGood = false
    index = 0
    orig_atoms = deepcopy(atoms)
    while !allGood
        index += 1
        println("E-max", E_max)
        a_rates = walk_single_walker!(atoms,model,NS.walker_params,E_max)
        display(a_rates)
        display(NS.walker_params)
        rates_good = []
        for key in keys(a_rates)
            rate = a_rates[key][2]/a_rates[key][1]
            push!(rates_good,adjust_step_sizes!(NS.walker_params,key,rate))

        end
        println("rates good")
        println(rates_good)
        allGood =  all(rates_good)
        ##println(allGood)
        atoms = deepcopy(orig_atoms)
        if index > 2000
            error("Taking too long")
        end
    end
        

#    possible = [do_cell_volume_step, do_cell_shear_step, do_cell_stretch_step, do_atoms_step]
#    # Loop over the number of random walks to take.
#    n_try = walk_params.n_single_walker_steps
#    n_accept = 0
#    idx = 0
#    acceptance_rates = Dict("stretch"=>0//n_try,"shear"=>0//n_try,"volume"=>0//n_try, "atoms"=>0//n_try)
#
#    orig_config = deepcopy(atoms)
##    volume_accept_rate = 0//n_try
##    while volume_accept_rate > 1//2  || volume_accept_rate < 1//4
##        n_accept = 0
##        for i = 1:walk_params.n_single_walker_steps
##            (__,n_accept_iter) = do_cell_volume_step(atoms,model,walk_params,E_max )
##            n_accept += n_accept_iter
##        end
##        volume_accept_rate = n_accept/n_try
##        if volume_accept_rate > 1//2
##            walk_params.volume_step_size *= 1.05
##        elseif volume_accept_rate < 1//4
##            walk_params.volume_step_size *= 0.95
##        end
##        println("Current step size volume: ", walk_params.volume_step_size)
##        println("Current accept rate: ", volume_accept_rate)
##    end
##
##
##    shear_accept_rate = 0//n_try
##    while shear_accept_rate > 1//2  || shear_accept_rate < 1//4
##        n_accept = 0
##        for i = 1:walk_params.n_single_walker_steps
##            (__,n_accept_iter) = do_cell_shear_step(atoms,model,walk_params,E_max )
##            n_accept += n_accept_iter
##        end
##        shear_accept_rate = n_accept/n_try
##        if shear_accept_rate > 1//2
##            walk_params.shear_step_size *= 1.05
##        elseif shear_accept_rate < 1//4
##            walk_params.shear_step_size *= 0.95
##        end
##        println("Current step size shear: ", walk_params.shear_step_size)
##        println("Current accept rate: ", shear_accept_rate)
##    end
##
##    stretch_accept_rate = 0//n_try
##    while stretch_accept_rate > 1//2  || stretch_accept_rate < 1//4
##        n_accept = 0
##        for i = 1:walk_params.n_single_walker_steps
##            (__,n_accept_iter) = do_cell_stretch_step(atoms,model,walk_params,E_max )
##            n_accept += n_accept_iter
##        end
##        stretch_accept_rate = n_accept/n_try
##        if stretch_accept_rate > 1//2
##            walk_params.stretch_step_size *= 1.05
##        elseif stretch_accept_rate < 1//4
##            walk_params.stretch_step_size *= 0.95
##        end
##        println("Current step size stretch: ", walk_params.stretch_step_size)
##        println("Current accept rate: ", stretch_accept_rate)
##    end    
#
#
#
#    atom_accept_rate = 0//n_try
#    index = 0
#    while atom_accept_rate > 1//2  || atom_accept_rate < 1//4
#        index += 1
#        n_accept = 0
#        for i = 1:walk_params.n_single_walker_steps
#            (__,n_accept_iter) = do_atoms_step(atoms,model,walk_params,E_max, KE_max = walk_params.KE_max )
#            n_accept += n_accept_iter
#        end
#        atom_accept_rate = n_accept//n_try
#        if atom_accept_rate > 1//2
#            walk_params.MC_atom_step_size *= 1.05
#        elseif atom_accept_rate < 1//4
#            walk_params.MC_atom_step_size *= 0.95
#        end
#        atoms = deepcopy(orig_config)
#        println("Current step size atom: ", walk_params.MD_time_step)
#        println("Current accept rate: ", atom_accept_rate, "<--------------------------------------------------------->")
#        if index > 100
#            error("Too long")
#        end
#    end    
            #    for iWalk in 1:walk_params.n_single_walker_steps

   # for i =1:10 # How many steps should we take to do this?
   #     idx += 1
#  #      println("Current Energy", ase.eval_energy(atoms,model))
   #     move = rand(1:length(possible))
   #     (tried,accepted) = possible[move](atoms,model,walk_params,E_max)
   #     acceptance_rates[string(possible_move)[9:end-5]] += accepted//n_try
#  #      n_try += n_iter_try
#  #      n_accept += n_iter_accept
   #     println("Walk Type: ", move)
   #     println("Index: ", idx)
   #     println("n_accept: ", n_accept)
   #     println("n_try: ", n_try)
   #     if idx > 300
   #         error("Too long.. Stopping")
   #     end
   #     
   #     
   # end

end

function reverse_sort_energies(NS)
    energies = zeros(MVector{length(NS.walkers),Float64})
    perms = zeros(MVector{length(NS.walkers),Int})
    for (iwalker, walker) in enumerate(NS.walkers)
        energies[iwalker] = walker.energies[2]
    end

    sortperm!(perms,energies)
    return perms
#    return reverse(perms)
        

end

function run_NS(NS::NS,LJ::LennardJones.model)

    V = (NS.n_walkers - NS.n_cull + 1)/(NS.n_walkers + 1)
    cDir = pwd()

    i = 1
    io = open(joinpath(cDir,"NS.out"),"w")
    write(io, "N_walkers = " * string(NS.n_walkers) * "\n")
    write(io, "N_cull = " * string(NS.n_cull) * "\n")
    write(io, "N_steps_per_walker = " * string(NS.n_iter) * "\n")
    write(io, "eps = " * string(NS.eps) * "\n")
    perms = zeros(MVector{length(NS.walkers),Int})
#    while V > NS.eps
    for i= 1:100
        println(i)
        println("V = ", V)
        ## Find the top Kr highest energies
        
        sortperm!(perms,[i.energies[2] for i in NS.walkers])  # 4 allocations
        perms = reverse(perms) # 1 allocation
        #perms =  reverse(sortperm([i.energies[2] for i in NS.walkers]))
       # println(sort([i.energyPerAtomModel for i in NS.walkers]))
        E_max = NS.walkers[perms[NS.n_cull]].energies[2]  # 3 allocations
        # Which configs need to be thrown out.
        forDelete = perms[1:NS.n_cull]  # 2 allocations
        # Which configs can be kept
        keeps = perms[NS.n_cull + 1: end]  # 2 allocations
        println("Energy cutoff")
        display(E_max)
        println("PE: ", ase.eval_energy(NS.walkers[perms[NS.n_cull]],LJ,do_KE = false))
        println("KE: ", ase.eval_KE(NS.walkers[perms[NS.n_cull]]))
        println("Total: ", ase.eval_energy(NS.walkers[perms[NS.n_cull]],LJ))

        if i %12 == 0  # 12 is pretty arbitrary.. Need a better way to see if need to re-tune
            println("Stopping to retune step sizes")
            energies_before = [ase.eval_energy(walker, LJ, force_recalc = true) for walker in NS.walkers]
            tune_step_sizes!(NS,LJ)
            energies_after = [ase.eval_energy(walker, LJ, force_recalc = true) for walker in NS.walkers]
            if !isapprox(energies_before,energies_after, atol=1e-5)
                error("tune step sizes changed configurations")
            end
            println("Done with tune up------------------------------------------------->")
        end
#        display(E_max - ase.eval_KE(NS.walkers[perms[NS.n_cull]]))
        for replace_walker in forDelete
            #Copy one of the configs that didn't get thrown out as the starting point
            #println("Initializing random walk to replace configuration ", i)
            #display(sample(keeps))
            NS.walkers[replace_walker] = deepcopy(NS.walkers[sample(keeps,1)[1]])  # ~ 50 allocations
            #println("E-max for this walker: ", E_max)
            #println("Starting energy of this walker: (should be lower than E-max)", ase.eval_energy(NS.walkers[replace_walker],LJ,force_recalc = true))
            walk_single_walker!(NS.walkers[replace_walker],LJ,NS.walker_params,E_max)
            #return nothing            
            #println("Ending energy of this walker: (should be lower than E-max)", ase.eval_energy(NS.walkers[replace_walker],LJ,force_recalc = true))
#            randomWalk!(NS.walkers[replace_walker],LJ,energyCutoff,NS.n_iter)
        end
        i += 1
        V = ((NS.n_walkers - NS.n_cull + 1)/(NS.n_walkers + 1))^i
        write(io,string(V) * " ")
        write(io,string(E_max) * " \n")

    end

    close(io)
end


function walk_single_walker!(atoms::ase.atoms, model, walk_params::NS_walker_params,E_max)
    possible = (do_atoms_step,do_cell_shear_step,do_cell_stretch_step,do_cell_volume_step)#,do_atoms_step)#, do_cell_stretch_step)#, do_atoms_step]
    # Loop over the number of random walks to take.
    n_try = 0
    n_accept = 0
    idx = 0
#    stretch_acceptance = ()
    atomsCount = 1
    typeLookup = Dict( 1 => "atoms", 2=>"shear",3=>"stretch",4=>"volume")
    acceptance_rates = Dict{String,Tuple{Int64,Int64}}("stretch"=>(0,0),"shear"=>(0,0),"volume"=>(0,0), "atoms"=>(0,0))
    for iWalk in 1:walk_params.n_single_walker_steps
        idx += 1
        move = rand(1:length(possible))
        #if move == 1
        #    atomsCount += 1
        #end
#        println("Doing move: ", typeLookup[move])
#        before_energy = ase.eval_energy(atoms,model,force_recalc = true)
#
#        if before_energy > E_max
#            println("Energy before move is greater than E_max")
#            println("iWalk = $iWalk")
#            println(before_energy)
#            println(E_max)
#            error("This should not have happened!")
#        end
        if move == 1
#            continue
            tried,accepted = do_atoms_step(atoms,model,walk_params,E_max)
        elseif move == 2
            tried,accepted = do_cell_shear_step(atoms,model,walk_params,E_max)
        elseif move == 3
            tried,accepted = do_cell_stretch_step(atoms,model,walk_params,E_max)
        elseif move == 4
            tried,accepted = do_cell_volume_step(atoms,model,walk_params,E_max)
        else
            error("out of bounds")
        end
            #        tried,accepted = possible[move](atoms,model,walk_params,E_max)
        #if idx == 200
        #    println("atomscount = $atomsCount")
        #    return nothing
        #end
#        after_energy = ase.eval_energy(atoms,model,force_recalc = true)
#        if after_energy > E_max
#            println("iWalk = $iWalk")
#            println(after_energy)
#            println(E_max)
#            error("This should not have happened!!!")
#        end
        acceptance_rates[typeLookup[move]] = (acceptance_rates[typeLookup[move]][1] + tried, acceptance_rates[typeLookup[move]][2] + accepted)
        if idx > 10000
            error("Too long.. Stopping")
        end
    end
    return acceptance_rates
end

function do_atoms_step(atoms::ase.atoms, model,walk_params,E_max) #atoms::ase.atoms, model,walk_params,E_max)


    if walk_params.atom_algorithm == "GMC"
        return ase.do_GMC!(atoms,model,walk_params,E_max)
    elseif walk_params.atom_algorithm == "MC"
        return ase.do_MC!(atoms,model,walk_params,E_max)
    elseif walk_params.atom_algorithm == "MD"
        return ase.do_MD!(atoms,model,walk_params,E_cutoff = E_max,KE_cutoff = walk_params.KE_max)
    end
    
end

function do_cell_shear_step(atoms::ase.atoms, model,walk_params,E_max;execute = true,check_energy = true)
  #  println("Doing cell shear")
    (p,T) = propose_shear_step(atoms,walk_params.shear_step_size)
    if do_cell_step!(atoms,p,T,E_max,walk_params,model,execute = execute,check_energy=check_energy)
        return (1,1)
    else
        return (1,0)
    end
end

function do_cell_stretch_step(atoms::ase.atoms, model,walk_params,E_max;execute=true,check_energy = true)
#    println("Doing cell stretch")
    p,T = propose_stretch_step(atoms,walk_params.stretch_step_size)
#    p = 1.0
#   T = SMatrix{3,3,Float64,9}(UniformScaling(1.0))
#    return T
    if do_cell_step!(atoms,p,T,E_max,walk_params,model,execute = execute,check_energy=check_energy)
        return 1,1
    else
        return 1,0
    end
end


function do_cell_volume_step(atoms::ase.atoms, model,walk_params,E_max;execute=true,check_energy = true)
#    println("Doing cell volume")
    (p,T) = propose_volume_step(atoms,walk_params.volume_step_size)
    
    if do_cell_step!(atoms,p,T,E_max,walk_params,model,execute = execute,check_energy=check_energy)
        return 1,1
    else
        return 1,0
    end
end

function do_cell_step!(atoms,p_accept,T,E_max,ns_params,model;execute = true,check_energy = true)
    #nAtoms = length(atoms.positions) # Don't need this...
    if (p_accept < 1.0) && (rand() > p_accept)
#        println("rejected for probability")
        return false   
    end

    new_cell = T * SMatrix(atoms.lVecs)
#    println("new cell right after definition")
#    display(T)
#    display(new_cell)

#    new_cell = SMatrix{3,3,Float64,9}(undef)  # Allocation
#    mul!(new_cell,atoms.lVecs,T)
#    newCell = atoms.lVecs * T
    newVolume = ase.cell_volume(new_cell * atoms.latpar)#atoms.latpar .* new_cell)  # Allocation from the multiplication ... fixed with mul!
#    println("Proposed Volume: ", newVolume)
#    println("Current Volume", ase.cell_volume(atoms))
#    println("Volume limit: ",ns_params.max_volume_per_atom *nAtoms)
    if newVolume > ns_params.max_volume_per_atom * atoms.nAtoms
#        println("Volume too large")
#        println("new volume: ", newVolume)
#        println("old volume: ", ase.cell_volume(atoms.latpar * atoms.lVecs))
#        println("Max Volume: ", ns_params.max_volume_per_atom * nAtoms)
        return false
    end
    if ase.min_aspect_ratio(new_cell) < ns_params.min_aspect_ratio

 #       println("Aspect ratio too small")
 #       println("new aspect ratio: ", ase.min_aspect_ratio(new_cell))
 #       println("old aspect ratio: ", ase.min_aspect_ratio(atoms.lVecs))
 #       println("Max aspect ratio: ", ns_params.min_aspect_ratio)
        return false
    end
    if check_energy
#        println("checking energy")
#         oldCell = deepcopy(atoms.lVecs)
#        oldPositions = atoms.positions  # I don't think I need this. set_cell! should take care of the undoing of positions
   #     println("before")
    #    display(ase.eval_energy(atoms,model))
    
#        atoms.lVecs = new_cell
    #    println("before set cell")
    #    display(atoms.lVecs)
  #      ase.DirectToCartesian!(atoms)
  #      display(atoms.positions)
#        println("setting it to this...")
#        display(new_cell)
        ase.set_cell!(atoms,new_cell,scale_atoms = true)
  #      println("after set cell")
#        display(atoms.lVecs)
  #      display(atoms.positions)
 #       println("after ")
 #       display(ase.eval_energy(atoms,model,force_recalc = true))
        newEnergy = ase.eval_energy(atoms,model,force_recalc = true)
  #      println("New Energy right before energy check: ", newEnergy)
        if newEnergy < E_max
 #           println("Accepted")
            # For tuning the step sizes we want to check to see if the move will be accepted 
            # but not actually do it.
            if execute
                atoms.energies[2] = newEnergy
            else
                old_cell = inv(T) * new_cell
                ase.set_cell!(atoms,old_cell,scale_atoms = true)
            end
            return true
        else
 #           println("Rejected for energy")
     #       println("New Energy:", newEnergy)
     #       println("E-max: ", E_max)
            old_cell = inv(T) * new_cell
            ase.set_cell!(atoms,old_cell,scale_atoms = true)
#            println("after resettng to original cell")
#            display(atoms.lVecs)
            #            atoms.positions .= oldPositions # Can I just set cell with scale_atoms = true again?
            return false
        end
    else
        # Not checking for energy, passed the volume and aspect ratio test, so it's a good step
        if execute
            ase.set_cell!(atoms,new_cell,scale_atoms = true)
        end
        return true
    end

end


function propose_volume_step(atoms::ase.atoms, stepsize::Float64)
    dV = rand(Normal(0,stepsize * atoms.nAtoms))
    # Get the current volume
    orig_V = ase.cell_volume(atoms)
    # Calculate the new volume
    new_V = orig_V + dV
    if new_V < 0
        println("Negative volume encountered, resetting to positive value.")
        new_V = abs(new_V)
    end
    # Calculate the scale factor

    transform = SMatrix{3,3,Float64,9}(UniformScaling((new_V / orig_V)^(1.0/3.0)))
    
    # This skews the probability towards larger volumes..
    # Because if new_V > orig_V, then the probability of accepting the move is 1.0, if new_V < orig_V, then the probability of accepting the move is less than 1.0. (Note... I'm still not sure why we want to do this. Waiting on a response from Livia)
    p_accept = minimum((1.0, (new_V/orig_V)^atoms.nAtoms))
    
    return p_accept, transform
end


function propose_shear_step(atoms::ase.atoms,stepsize::Float64)
    rmv_vec = rand(1:3)
    remaining_vecs = rmv_vec == 1 ? (2,3) : rmv_vec == 2 ? (1,3) : (1,2)
    #new_cell = copy(atoms.lVecs)
   # new_cell = MMatrix{3,3,Float64,9}(undef)
   # new_cell .= atoms.lVecs
#    new_cell = buffer
#    new_cell .= atoms.lVecs
    v1 = atoms.lVecs[:,remaining_vecs[1]]
    v2 = atoms.lVecs[:,remaining_vecs[2]]
    norm_v1 = sqrt(dot(v1,v1))
    v1 = v1 / norm_v1
    v2 = v2 -  v1 * dot(v1, v2)  # Make v2 orthogonal to v1
    norm_v2 = sqrt(dot(v2,v2))
    v2 = v2/norm_v2
    if abs(dot(v1, v2)) > 1e-4
        println("v1 and v2 are not orthogonal, something is wrong.")
    end
#    rv = MVector{2,Float64}(undef)
#    rand!(Xoshiro(), Normal(0, stepsize), rv)
#    rv1, rv2 = rv
    rv1 = rand(Normal(0,stepsize))
    rv2 = rand(Normal(0,stepsize))
    new_vec = atoms.lVecs[:,rmv_vec] + rv1 * v1 + rv2 * v2
    l1 = atoms.lVecs[:,1]
    l2 = atoms.lVecs[:,2]
    l3 = atoms.lVecs[:,3]
    if rmv_vec == 1
        new_cell = hcat(new_vec,l2,l3)
    elseif rmv_vec == 2
        new_cell = hcat(l1,new_vec,l3)
    elseif rmv_vec == 3
        new_cell = hcat(l1,l2,new_vec)
    else
        error("BAD")
    end
    #rv1,rv2 = rand(Normal(0,stepsize),2)
#    new_cell[:,rmv_vec] .= atoms.lVecs[:,rmv_vec] .+ rv1 .* v1 .+ rv2 .* v2
#    println("after", atoms.lVecs)   
    inv_lVecs = inv(SMatrix{3,3,Float64,9}(atoms.lVecs))
#    println("Here")
#    println(typeof(new_cell),typeof(inv_lVecs))
   #println("new cell shear")
   # display(new_cell)
   # println("before shear")
   # display(atoms.lVecs)
    transform = new_cell * inv_lVecs
   # println("transform right after definition")
   # display(transform)
    return 1.0,SMatrix{3,3,Float64,9}(transform)

end


function propose_stretch_step(atoms::ase.atoms,stepsize::Float64) 
    # Yes, I know that I don't need an atoms object to do this, but to stay consistent with the shear step, I'm going to keep it.


#    a = stepsize .* (2 .* rand(2) .- 1)  # Get 3 random floats
#    A = [a[1] 0  0
#         0  a[2] 0
#         0  0  1/a[1]/a[2]]
#    return A

    rmv_vec_one = rand(1:3)
    rmv_vec_two = rand(1:3)
    if rmv_vec_one == rmv_vec_two
        #println("found duplicate vectors, regenerating rmv_vec_two")
        rmv_vec_two = (rmv_vec_two + 1) % 3 + 1
    end
    remaining_vec = ((5 - (rmv_vec_one + rmv_vec_two)) % 3) + 1 #[x for x in 1:3 if x != rmv_vec_one && x != rmv_vec_two][1]

    rv = rand(Normal(0,stepsize))
    if remaining_vec == 1
        v1 = SVector{3,Float64}(1,0,0)
        v2 = exp(rv) * SVector{3,Float64}(0,1,0)
        v3 = exp(-rv) * SVector{3,Float64}(0,0,1)
    elseif remaining_vec == 2
        v1 = exp(rv) * SVector{3,Float64}(1,0,0)
        v2 = SVector{3,Float64}(0,1,0)
        v3 = exp(-rv) * SVector{3,Float64}(0,0,1)
    elseif remaining_vec == 3
        v1 = exp(rv) * SVector{3,Float64}(1,0,0)
        v2 = exp(-rv) * SVector{3,Float64}(0,1,0)
        v3 = exp(rv) * SVector{3,Float64}(0,0,1)
    else
        error("BAD")
    end
    transform = [v1 v2 v3 ]
#    transform = zeros(MMatrix{3,3,Float64,9})
#    transform[remaining_vec,remaining_vec] = 1.0
#    transform[rmv_vec_one,rmv_vec_one] = exp(rv)
#    transform[rmv_vec_two,rmv_vec_two] = exp(-rv)
    return 1.0, transform
end



function do_cell_shape_walk!(atoms::ase.atoms, walker_params)
    possibilities =[do_cell_shear_step,do_cell_stretch_step, do_cell_volume_step]
    # I don't need an energy model to walk the cell shape, but the routine I want to use to do it
    # requires it, so I'm just initializing an empty model.
    empty_model = LennardJones.initialize_empty_model()
    #n_try = 0
    #n_accept = 0
    n_try = walker_params.n_single_walker_steps
    acceptance_rates = Dict("stretch"=>0//n_try,"shear"=>0//n_try,"volume"=>0//n_try)
    for i =1:walker_params.n_single_walker_steps
        shuffle!(possibilities)
        for possible_move in possibilities
            #println(walker_params)
            
            (tried,accepted) = possible_move(atoms,empty_model,walker_params,1000.0,check_energy = false)
            acceptance_rates[string(possible_move)[9:end-5]] += accepted//n_try
#            n_accept += accepted
#            n_try += tried
            #            if (p_accept < 1.0) & (rand() > p_accept)
#                atoms.lVecs = transform * atoms.lVecs
#                ase.DirectToCartesian!(atoms)
           #     println("Accepted")
#            end
        end
    end
    return atoms
#    println("tried: ", walker_params.n_single_walker_steps, "Accepted: ",acceptance_rates)
end








end