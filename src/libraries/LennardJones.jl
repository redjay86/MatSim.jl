#module LJ


function forceOnSingleParticle(positions::Array{SVector{2,Float64},1},particle::SVector{2,Float64},boxSize::Float64):: SVector{2,Float64}
    fVec = SVector(0,0)
    #modifiedPos = SVector{2,Float64}
    #diffVec = SVector{2,Float64}
    for i=1:size(positions,1)
        diffVec = particle - positions[i]

        if abs(diffVec[1]) > boxSize/2 && abs(diffVec[2]) > boxSize/2
            modifiedPos = positions[i] + boxSize * SVector(sign(diffVec[1]) , sign(diffVec[2]))
        elseif abs(diffVec[1]) > boxSize/2
            modifiedPos = positions[i] + boxSize * SVector(sign(diffVec[1]) , 0)
        elseif abs(diffVec[2]) > boxSize/2
            modifiedPos = positions[i] + boxSize * SVector(0 , sign(diffVec[2]))
        else
            modifiedPos = copy(positions[i])
        end
        diffVec = particle - modifiedPos
        distance = norm(diffVec)
        
        if distance > 0.5
            fVec += 24 * (2/distance^13 - 1/distance^7) * diffVec/distance
        end
    end
    return fVec

end




function singleAtomEnergy(crystal::Crystal,centerAtom::Vector{Float64}, centerType:: Integer, loopBounds::Vector{Int64},params::Matrix{Float64},cutoff:: Float64)
    ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    for (iNeighbor,aType) in enumerate(crystal.atomicBasis), neighboratom in aType  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + [i, j, k]
            newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < cutoff && !isapprox(r,0,atol = 1e-3)
                #If the neighbor atom is inside the unit cell, then its going to be
                # double counted at some point when we center on the other atom.  
                # So we count it as half each time.
                if all([i,j,k] .== 0 ) 
                    ljvals[iNeighbor + centerType - 1,1] +=  4 * params[iNeighbor + centerType - 1,1] * 1/2 * params[iNeighbor + centerType - 1,2]^6/r^6
                    ljvals[iNeighbor + centerType - 1,2] +=  4 * params[iNeighbor + centerType - 1,1] * 1/2 * params[iNeighbor + centerType - 1,2]^12/r^12
                else 
                    ljvals[iNeighbor + centerType - 1,1] += 4 * params[iNeighbor + centerType - 1,1] * params[iNeighbor + centerType - 1,2]^6/r^6
                    ljvals[iNeighbor + centerType - 1,2] += 4 * params[iNeighbor + centerType - 1,1] * params[iNeighbor + centerType - 1,2]^12/r^12
                end
            end
        end
    end
    return ljvals
end

function totalEnergy(crystal::Crystal,cutoff::Float64;params = Float64[1. 1.; 1. 1.; 1. 1.])
    CartesianToDirect!(crystal)
    ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    
    loopBounds = convert.(Int64,cld.(cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] ))
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis), centerAtom in centerAtomType 
        centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
        ljvals .+= singleAtomEnergy(crystal,centerAtomC,iCenter,loopBounds,params,cutoff)    # Find the contribution to the LJ energy for this centering atom.
    end
    return ljvals
end

function totalEnergy(crystal::Crystal,params::Array{Float64,2})
    energy = sum(-params[:,1] .* params[:,2] .^6 .* crystal.ljvals[:,1] + params[:,1] .* params[:,2] .^12 .* crystal.ljvals[:,2] )
    return energy
end


function loglik(data::DataSet,μ,σ)
    n = length(data.crystals)
    thesum =  - 1/(2 * σ^2) * sum([(i.energyFP - totalEnergy(i,μ))^2   for i in data.crystals])
    return thesum

end


#end