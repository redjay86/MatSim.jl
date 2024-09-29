
function initializeSimulation(inputs,species,model)
    parsedlVecs = [parse.(Float64,split(x[1])) for x in inputs["lVecs"]]
    lVecs = hcat(parsedlVecs...)
    configs = [buildRandom(inputs["lPar"],lVecs,inputs["nAtoms"],inputs["minSep"],species) for i in 1:inputs["K"]]
    for i in configs
        i.modelEnergy = MatSim.totalEnergy(i,model)
    end
    return NS(inputs["K"],inputs["Kr"],inputs["L"],inputs["eps"],inputs["walkMethod"],configs)
end

function nestedSampling(NS::NS,LJ::LJ)

    V = (NS.K - NS.Kr + 1)/(NS.K + 1)

    i = 1
    while V > NS.eps
        println(i)
        ## Find the top Kr highest energies
        sEnergies =  reverse(sortperm([i.modelEnergy for i in NS.configs]))
        energyCutoff = NS.configs[sEnergies[NS.Kr]].modelEnergy
        # Which configs need to be thrown out.
        forDelete = sEnergies[1:NS.Kr]
        # Which configs can be kept
        keeps = sEnergies[NS.Kr + 1: end]

        println("Energy cutoff")
        display(energyCutoff)
        for i in forDelete
            #Copy one of the configs that didn't get thrown out as the starting point
            println("Initializing random walk to replace configuration ", i)
            display(sample(keeps))
            NS.configs[i] = NS.configs[sample(keeps,1)[1]]
            randomWalk!(NS.configs[i],LJ,energyCutoff,NS.L)
        end
        i += 1
        V = ((NS.K - NS.Kr + 1)/(NS.K + 1))^i

    end

end

# Perform a random walk on a single configuration subject to the constraint that the total energy be less than the cutoff
# Not walking using the force vector (GMC), just random walks.  This may not work well for some systems.
function randomWalk!(config::Crystal,model, energyCutoff::Float64, nWalk::Int64)
    # Loop over the number of random walks to take.
    for iWalk in 1:nWalk
        @printf "Step in random walk %3i\n" iWalk
        #Loop over all of the atoms in the simulation.
        for (iType,atomType) in enumerate(config.atomicBasis), (iAtom,atom) in enumerate(atomType)
#            @printf "Atom being moved. Type: %3i Number: %3i\n" iType iAtom 
            # Get a random displacement vector
            randDisplace = (rand(3).-0.5)*0.1
            config.atomicBasis[iType][iAtom] .+= randDisplace
            newEnergy = totalEnergy(config,model)
            # If the move resulted in a higher energy, undo the move and go to the next atom.
            if newEnergy > energyCutoff
 #               @printf "Rejected\n"
                config.atomicBasis[iType][iAtom] .-= randDisplace
            else # Otherwise, update the energy 
#                @printf "Accepted\n"
                config.modelEnergy = newEnergy
            end     
        end   
#        randDisplacement = [[ (rand(3).-0.5)*0.1 for i =1:config.nType[j]] for j = 1:config.order]
#        config.atomicBasis .+= randDisplacement  # Use the displacement to move this atom.

       # println(totalEnergy(config,model))
    end
end

# Galilean Monte Carlo
 function GMC(config::Crystal,nWalk::Int64)

    oldEnergy = sum(totalEnergy(config,LJcutoff,params = LJparams))
    newEnergy = 1e6
    dt = 0.1
    
    velocities = [[ (rand(3).-0.5)*0.1 for i =1:config.aType[j]] for j = 1:config.order]
    for iWalk = 1:nWalk

        while newEnergy > oldEnergy
            
 #           config.atomicBasis += velocities * dt
 #           if any()

        end
    end

end



#nType = 3
#aType = [5, 2 , 3]
#test = [[(rand(3).-0.5)*0.1 for i =1:aType[j]] for j = 1:nType]
#display(test[2])
#
#M = 8
#m = 1.2
#d = 1.0
#d^2/2 *(1/2 * M + m)