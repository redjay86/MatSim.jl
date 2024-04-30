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




function singleAtomEnergy(LJ::LJ,crystal::Crystal,centerAtom::Vector{Float64}, centerType:: Integer, loopBounds::Vector{Int64})
    #ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    totalEnergy = 0
    for (iNeighbor,aType) in enumerate(crystal.atomicBasis), neighboratom in aType  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + [i, j, k]
            newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < LJ.cutoff && !isapprox(r,0,atol = 1e-3)
                #If the neighbor atom is inside the unit cell, then its going to be
                # double counted at some point when we center on the other atom.  
                # So we count it as half each time.
                if all([i,j,k] .== 0 ) 
                    totalEnergy +=  4 * LJ.params[iNeighbor + centerType - 1,1] * 1/2 * LJ.params[iNeighbor + centerType - 1,2]^6/r^6
                    totalEnergy +=  4 * LJ.params[iNeighbor + centerType - 1,1] * 1/2 * LJ.params[iNeighbor + centerType - 1,2]^12/r^12
                else 
                    totalEnergy += 4 * LJ.params[iNeighbor + centerType - 1,1] * LJ.params[iNeighbor + centerType - 1,2]^6/r^6
                    totalEnergy += 4 * LJ.params[iNeighbor + centerType - 1,1] * LJ.params[iNeighbor + centerType - 1,2]^12/r^12
                end
            end
        end
    end
    return totalEnergy
end

function singleAtomDistances!(distMat::Matrix{Float64},crystal::Crystal,LJ::LJ,centerAtom::Vector{Float64}, centerType:: Integer, loopBounds::Vector{Int64})
#    ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    for (iNeighbor,aType) in enumerate(crystal.atomicBasis), neighboratom in aType  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + [i, j, k]
            newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < LJ.cutoff && !isapprox(r,0,atol = 1e-3)
                #If the neighbor atom is inside the unit cell, then its going to be
                # double counted at some point when we center on the other atom.  
                # So we count it as half each time.
                if all([i,j,k] .== 0 ) 
                    distMat[iNeighbor + centerType - 1,1] +=  1/2 * 1/r^6
                    distMat[iNeighbor + centerType - 1,2] +=  1/2 * 1/r^12
                else 
                    distMat[iNeighbor + centerType - 1,1] += 1/r^6
                    distMat[iNeighbor + centerType - 1,2] += 1/r^12
                end
            end
        end
    end
    return distMat
end

function totalDistances(crystal::Crystal,LJ::LJ)
    CartesianToDirect!(crystal)
    distMat = zeros(size(LJ.params)...)#zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    
    loopBounds = convert.(Int64,cld.(LJ.cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] ))
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis), centerAtom in centerAtomType 
        centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
        singleAtomDistances!(distMat,crystal,LJ,centerAtomC,iCenter,loopBounds)    # Find the contribution to the LJ energy for this centering atom.
    end
    return distMat
end

function totalEnergy(crystal::Crystal,LJ::LJ)
    if !all(crystal.ljvals .== 0.0)
        totalEnergy = sum(-LJ.params[:,1] .* LJ.params[:,2] .^6 .* crystal.ljvals[:,1] + LJ.params[:,1] .* LJ.params[:,2] .^12 .* crystal.ljvals[:,2] ) 
        return totalEnergy
    end
    println("Doing it the hard way!")
    CartesianToDirect!(crystal)
    totalEnergy = 0 
    loopBounds = convert.(Int64,cld.(LJ.cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] ))
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis), centerAtom in centerAtomType 
        centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
        totalEnergy += singleAtomEnergy(LJ,crystal,centerAtomC,iCenter,loopBounds)    # Find the contribution to the LJ energy for this centering atom.
    end
    return totalEnergy
end

#function totalEnergy(crystal::Crystal,params::Array{Float64,2})
#    energy = sum(-params[:,1] .* params[:,2] .^6 .* crystal.ljvals[:,1] + params[:,1] .* params[:,2] .^12 .* crystal.ljvals[:,2] )
#    return energy
#end


function loglik(data::DataSet,energyModel,σ)
    n = length(data.crystals)
    thesum =  - 1/(2 * σ^2) * sum([(i.energyFP - totalEnergy(i,energyModel))^2   for i in data.crystals])
    return thesum

end

function initializeLJ(settings::Dict)

    order = settings["order"]
    cutoff = settings["cutoff"]
    nInteractions = Int(order*(order + 1)/2)
    params = ones(nInteractions,2)
    return LJ(order,params,cutoff)
end

#end