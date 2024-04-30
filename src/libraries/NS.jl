
function initializeSimulation(inputs,species,model)
    parsedlVecs = [parse.(Float64,split(x[1])) for x in inputs["lVecs"]]
    lVecs = hcat(parsedlVecs...)
    configs = [buildRandom(inputs["lPar"],lVecs,inputs["nAtoms"],inputs["minSep"],species) for i in 1:inputs["K"]]
    for i in configs
        i.energyPred = totalEnergy(i,model)
    end
    return NS(inputs["K"],inputs["Kr"],inputs["L"],inputs["eps"],inputs["walkMethod"],configs)
end

function nestedSampling(NS::NS)

    V = (NS.K - NS.Kr + 1)/(NS.K + 1)

    i = 1
    while V > NS.eps
        ## Find the top Kr highest energies

        i += 1
        V = ((NS.K - NS.Kr + 1)/(NS.K + 1))^i

    end

end
# Perform a random walk on a single configuration subject to the constraint that the total energy be less than the cutoff
function randomWalk(config::Crystal,model, energyCutoff::Float64, nWalk::Int64)
    # Loop over the number of random walks to take.
    oldEnergy = totalEnergy(config,model)
    newEnergy = 1e6
    for iWalk in 1:nWalk
        #Loop over all of the atoms in the simulation.
        # Loop until you find a displacement that is downhill
        while newEnergy > oldEnergy
            randomDisplacment = [[ (rand(3).-0.5)*0.1 for i =1:config.aType[j]] for j = 1:config.order]
            config.atomicBasis .+= randDisplacement  # Use the displacement to move this atom.
            newEnergy = totalEnergy(config,model)
        end
        # Once I've found a config with a lower energy, set the new energy to be the old in preparation for the next step in the random walk.
        oldEnergy = newEnergy  
        newEnergy = 1e6
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