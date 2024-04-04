
function initializeSimulation(lPar,lVecs,minSep,nAtoms,nConfigs,nWalks,V_ϵ)
    configs = [buildRandom(lPar,lVecs,nAtoms,minSep) for i in 1:nConfigs]
    return NS(lPar,nConfigs,nWalks,V_ϵ,nAtoms,configs)
end

function randomWalk(config::Crystal,nWalk::Int64,LJcutoff::Float64,LJparams::Array{Float64,2})
    # Loop over the number of random walks to take.
    oldEnergy = sum(totalEnergy(config,LJcutoff,params = LJparams))
    newEnergy = 1e6
    for iWalk in 1:nWalk
        #Loop over all of the atoms in the simulation.
        # Loop until you find a displacement that is downhill
        while newEnergy > oldEnergy
            randomDisplacment = [[ (rand(3).-0.5)*0.1 for i =1:config.aType[j]] for j = 1:config.order]
            config.atomicBasis .+= randDisplacement  # Use the displacement to move this atom.
            newEnergy = sum(totalEnergy(config,LJcutoff,params = LJparams))
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