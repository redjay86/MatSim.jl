module nSampling

using Crystal
using LennardJones
using DataSets:DataSet
using StaticArrays
using LinearAlgebra
using Distributions

struct NS
    K:: Int
    Kr:: Int
    L:: Int
    eps:: Float64
    method:: String
    configs::Vector{Crystal.config}
end


function randomShearMatrix(bound)
    a = bound .* (2 .* rand(3) .- 1)  # Get 3 random floats
    A = [1. a[1] a[2]
         0   1.  a[3]
         0   0   1.0]
    return A

end


function randomStretchMatrix(bound)
    a = bound .* (2 .* rand(2) .- 1)  # Get 3 random floats
    A = [a[1] 0  0
         0  a[2] 0
         0  0  1/a[1]/a[2]]
    return A
end
function initialize(inputs,species,model)
    parsedlVecs = [parse.(Float64,split(x[1])) for x in inputs["lVecs"]]
    lVecs = hcat(parsedlVecs...)
    configs = [Crystal.buildRandom(inputs["lPar"],lVecs,inputs["nAtoms"],inputs["minSep"],species) for i in 1:inputs["K"]]
    for i in configs
        i.energyPerAtomModel = LennardJones.totalEnergy(i,model)
    end
    return NS(inputs["K"],inputs["Kr"],inputs["L"],inputs["eps"],inputs["walkMethod"],configs)
end

function simulate(NS::NS,LJ::LennardJones.model)

    V = (NS.K - NS.Kr + 1)/(NS.K + 1)

    i = 1
    while V > NS.eps
        println(i)
        println("V = ", V)
        ## Find the top Kr highest energies
        sEnergies =  reverse(sortperm([i.energyPerAtomModel for i in NS.configs]))
       # println(sort([i.energyPerAtomModel for i in NS.configs]))
        energyCutoff = NS.configs[sEnergies[NS.Kr]].energyPerAtomModel  # Find the 
        # Which configs need to be thrown out.
        forDelete = sEnergies[1:NS.Kr]
        # Which configs can be kept
        keeps = sEnergies[NS.Kr + 1: end]

        println("Energy cutoff")
        display(energyCutoff)
        for i in forDelete
            #Copy one of the configs that didn't get thrown out as the starting point
            #println("Initializing random walk to replace configuration ", i)
            #display(sample(keeps))
            NS.configs[i] = NS.configs[sample(keeps,1)[1]]
            randomWalk!(NS.configs[i],LJ,energyCutoff,NS.L)
        end
        i += 1
        V = ((NS.K - NS.Kr + 1)/(NS.K + 1))^i

    end

end

# Perform a random walk on a single configuration subject to the constraint that the total energy be less than the cutoff
# Not walking using the force vector (GMC), just random walks.  This may not work well for some systems.
function randomWalk!(config::Crystal.config,model, energyCutoff::Float64, nWalk::Int64)
    # Loop over the number of random walks to take.
    for iWalk in 1:nWalk
        #@printf "Step in random walk %3i\n" iWalk
        #Loop over all of the atoms in the simulation.
        for (iType,atomType) in enumerate(config.atomicBasis), (iAtom,atom) in enumerate(atomType)
            #@printf "Atom being moved. Type: %3i Number: %3i\n" iType iAtom 
            # Get a random displacement vector
            randDisplace = (rand(3).-0.5)*0.01
#            println("before")
#            display(config.atomicBasis[iType][iAtom])
            config.atomicBasis[iType][iAtom] = (config.atomicBasis[iType][iAtom] + randDisplace) .% 1
#            println("after")
#            display(config.atomicBasis[iType][iAtom])
#            println(config.coordSys)
            newEnergy = LennardJones.totalEnergy(config,model) # Want total energies for NS
            # If the move resulted in a higher energy, undo the move and go to the next atom.
            if newEnergy > energyCutoff
                #println("Rejected")
                config.atomicBasis[iType][iAtom] -= randDisplace
            else # Otherwise, update the energy 
                #println("Accepted")
                config.energyPerAtomModel = newEnergy
            end     
        end   
#        randDisplacement = [[ (rand(3).-0.5)*0.1 for i =1:config.nType[j]] for j = 1:config.order]
#        config.atomicBasis .+= randDisplacement  # Use the displacement to move this atom.

       # println(totalEnergy(config,model))
    end
end

# Galilean Monte Carlo
 function GMC(config::Crystal.config,nWalk::Int64,model)
    initialConfig = deepcopy(config)
    oldEnergy = LennardJones.totalEnergy(config,model)
    Crystal.DirectToCartesian!(config)
    println("energy at start")
    println(oldEnergy)
    newEnergy = 1e6
    dt = 0.1
    
    velocities = [[zeros(SVector{3}) for i = 1:length(config.atomicBasis[j])] for j = 1:config.order]
    for j = 1:config.order, i = 1:1:length(config.atomicBasis[j])
        randVel = convert(SVector{3},(  2 * rand(3) .- 1.0))
        velocities[j][i] =randVel/norm(randVel)
    end
    forces = [[zeros(SVector{3}) for x = 1:length(config.atomicBasis[n])] for n = 1:length(config.atomicBasis)]
    display(velocities)
    for iWalk = 1:nWalk
        println("iWalk: ", iWalk)
        println(config.coordSys)
        display(config.atomicBasis[1][1])
        config.atomicBasis += velocities * dt  # Propogate the postions forward
        println("after movement")
        println(config.coordSys)
        display(config.atomicBasis[1][1])
        newEnergy = LennardJones.totalEnergy(config,model)  # Calculate the new energies
        Crystal.DirectToCartesian!(config)
        #println("new locations: ")
        #display(config.atomicBasis)
        println("newEnergy: ", newEnergy)
        if newEnergy > oldEnergy  # If we went uphill, we need to try and re-direct the velocites in the direction of the net force.
            println("Redirecting....----------------------------------->")
            for (iType,aType) in enumerate(config.atomicBasis), (iAtom,atom) in enumerate(aType)  #Loop over the different atom types.
                forces[iType][iAtom] = LennardJones.gradientForce(model,config,@SVector[iType,iAtom],@SVector[2,2,2])
            end
            nHat = [x ./ norm.(x) for x in forces]
#            println("velocity")
#            display(velocities[1][1])
#            display(nHat[1][1])
            for j = 1:config.order, i = 1:length(config.atomicBasis[j])
                @inbounds velocities[j][i] -= (2 * (velocities[j][i]' * nHat[j][i])) * nHat[j][i]
            end
#            println("after deflection")
 #           display(velocities[1][1])
#            dots = [[velocities[j][i]' * nHat[j][i] * nHat[j][i] for i =1:length(config.atomicBasis[j])] for j = 1:config.order]
            #println("dots")
            #display(dots)
 #           velocities = velocities - 2 * dots
        end
    end

    if newEnergy > oldEnergy
        return initialConfig
    else
        return config
    end
end



end