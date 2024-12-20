

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

function gradientForce(LJ,crystal,atom,loopBounds; eps = 1e-5)
    fVec = zeros(3)
    # Find x-component of force
#    println("before")
#    display(crystal.atomicBasis[atom[1]][atom[2]])
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] -= SVector{3,Float64}(eps,0,0)
    energyOne = MatSim.totalEnergy(crystal,LJ)
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] += SVector{3,Float64}(2 * eps,0,0)
    energyTwo = MatSim.totalEnergy(crystal,LJ)
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] -= SVector{3,Float64}(eps,0,0)  # Put it back where it was.
    fVec[1] = -(energyTwo - energyOne)/(2 * eps)
#    println("energies after 1")
#    display(energyTwo)
#    display(energyOne)
    # Find y-component of force
    crystal.atomicBasis[atom[1]][atom[2]] -= SVector{3,Float64}(0,eps,0)
    energyOne = MatSim.totalEnergy(crystal,LJ)
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] += SVector{3,Float64}(0,2 * eps,0)
    energyTwo = MatSim.totalEnergy(crystal,LJ)
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] -= SVector{3,Float64}(0,eps,0)  # Put it back where it was.
    fVec[2] = -(energyTwo - energyOne)/(2 * eps)

 #   println("energies after 2")
 #   display(energyTwo)
 #   display(energyOne)
    # Find z-component of force
    crystal.atomicBasis[atom[1]][atom[2]] -= SVector{3,Float64}(0,0,eps)  # Move atom downward
    energyOne = MatSim.totalEnergy(crystal,LJ)        # Calculate energy
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] += SVector{3,Float64}(0,0,2 * eps)  # Move atom upward
    energyTwo = MatSim.totalEnergy(crystal,LJ)           # Calculate energy
    MatSim.DirectToCartesian!(crystal)
    crystal.atomicBasis[atom[1]][atom[2]] -= SVector{3,Float64}(0,0,eps)  # Put it back where it was.
    fVec[3] = -(energyTwo - energyOne)/(2 * eps)      # Calculate component of gradient
 #   println("energies after 3")
 #   display(energyTwo)
 #   display(energyOne)

    return fVec
end

function singleAtomForce(LJ::LJ,crystal::Crystal,centerAtom::SVector{2,Int64}, loopBounds::SVector{3,Int64})
    #ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to k-nary case.
    MatSim.CartesianToDirect!(crystal)

    addVec = zeros(3)
    indices = zeros(2)
    fVec = SVector(0,0,0)
    for (iNeighbor,aType) in enumerate(crystal.atomicBasis), neighboratom in aType  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            addVec .= (i,j,k) 
            newAtom = neighboratom + addVec  #Periodic image of this atom
            newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)  # Convert to cartesian coordinate system
            r = newCart - MatSim.DirectToCartesian(crystal.latpar * crystal.lVecs,crystal.atomicBasis[centerAtom[1]][centerAtom[2]]) 
            if norm(r) < LJ.cutoff && !isapprox(norm(r),0,atol = 1e-3)
                println("Adding to force")
                indices = iNeighbor < centerAtom[1] ? @SVector[iNeighbor,centerAtom[1]] : @SVector[centerAtom[1],iNeighbor]
                fVec -=    12. * 4. * LJ.ϵ[indices[1],indices[2]] * LJ.σ[indices[1],indices[2]]^12/norm(r)^13 * r/norm(r)
                fVec +=    6. * 4. * LJ.ϵ[indices[1],indices[2]] * LJ.σ[indices[1],indices[2]]^6/norm(r)^7 * r/norm(r)
            end
        end
    end
    return fVec
end


function singleAtomEnergy(LJ::LJ,crystal::Crystal,centerAtom::SVector{3,Float64}, centerType:: Integer, loopBounds::SVector{3,Int64})
    #ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    totalEnergy = 0
    addVec = zeros(3)
    indices = zeros(2)
    for (iNeighbor,aType) in enumerate(crystal.atomicBasis), neighboratom in aType  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            addVec .= (i,j,k) 
#            addVec[2] = Float64(j)
#            addVec[3] = Float64(k)
#            addVec[1],addVec[2],addVec[3] .= Float64(i), Float64(j), Float64(k)
            newAtom = neighboratom + addVec
            newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < LJ.cutoff && !isapprox(r,0,atol = 1e-3)
                #If the neighbor atom is inside the unit cell, then its going to be
                # double counted at some point when we center on the other atom.  
                # So we count it as half each time.
                indices = iNeighbor < centerType ? @SVector[iNeighbor,centerType] : @SVector[centerType,iNeighbor]
                if all(isapprox.(addVec,0.0) ) 
                    totalEnergy -=  4. * LJ.ϵ[indices[1],indices[2]] * 1/2 * LJ.σ[indices[1],indices[2]]^6/r^6
                    totalEnergy +=  4. * LJ.ϵ[indices[1],indices[2]] * 1/2 * LJ.σ[indices[1],indices[2]]^12/r^12
                else 
                    totalEnergy -= 4. * LJ.ϵ[indices[1],indices[2]] * LJ.σ[indices[1],indices[2]]^6/r^6
                    totalEnergy += 4. * LJ.ϵ[indices[1],indices[2]] * LJ.σ[indices[1],indices[2]]^12/r^12
                end
            end
        end
    end
    return totalEnergy
end

function singleAtomDistances!(crystal::Crystal,LJ::LJ,centerAtom::SVector{3,Float64}, centerType:: Integer, loopBounds::SVector{3,Int64})
#    ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    for (iNeighbor,aType) in enumerate(crystal.atomicBasis), neighboratom in aType  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + @SVector[i,j,k]
            newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < LJ.cutoff && !isapprox(r,0.0,atol = 1e-3)
                #If the neighbor atom is inside the unit cell, then its going to be
                # double counted at some point when we center on the other atom.  
                # So we count it as half each time.
                # The LJ parameters are stored in the upper triangular portion of a matrix
                # The bottom triangle is redundant.. interaction between a and b is equivalent
                # to interaction between b and a.  So I sort the indices here so that the bottom
                # triangle of the matrix never gets updated, only the upper right.
                #indices[1],indices[2] = iNeighbor,centerType
                
                indices = iNeighbor < centerType ? @SVector[iNeighbor,centerType] : @SVector[centerType,iNeighbor]
                if all(@SVector[i,j,k] .== 0 ) 
                    crystal.r6[indices[1],indices[2]] +=  4. * 1.0/2.0 * 1.0/r^6
                    crystal.r12[indices[1],indices[2]] +=  4. * 1.0/2.0 * 1.0/r^12
                else 
                    crystal.r6[indices[1],indices[2]] += 4. * 1.0/r^6
                    crystal.r12[indices[1],indices[2]] += 4. * 1.0/r^12
                end
            end
        end
    end
    #return distMat
end

function totalDistances!(crystal::Crystal,LJ::LJ)
    CartesianToDirect!(crystal)
#    r6 = zeros(crystal.order,crystal.order)
#    r12 = zeros(crystal.order,crystal.order)
    
    eVals = eigvals(transpose(crystal.latpar .* crystal.lVecs) * (crystal.latpar .* crystal.lVecs))
    maxN = LJ.cutoff/sqrt(minimum(eVals))
    
    loopBounds = SVector{3,Int64}(convert.(Int64,cld.(LJ.cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)) )))
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis), centerAtom in centerAtomType 
        centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
        singleAtomDistances!(crystal,LJ,centerAtomC,iCenter,loopBounds)    # Find the contribution to the LJ energy for this centering atom.
        
    end
end

function totalEnergy(crystal::Crystal,LJ::LJ)
    if !all(crystal.r6 .== 0.0)
#        println("Doing it the easy way")
        totalEnergy = 0.0
        for i in eachindex(LJ.ϵ)
            totalEnergy += -LJ.ϵ[i] * LJ.σ[i]^6 * crystal.r6[i] + LJ.ϵ[i] * LJ.σ[i]^12 * crystal.r12[i]
        end
    else
 #       println("Doing it the hard way")
        CartesianToDirect!(crystal)
        totalEnergy = 0.0 
        loopBounds = SVector{3,Int64}(convert.(Int64,cld.(LJ.cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] )))
        # The outer two loops are to loop over different centering atoms.
        for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis), centerAtom in centerAtomType 
            centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
            totalEnergy += singleAtomEnergy(LJ,crystal,centerAtomC,iCenter,loopBounds)    # Find the contribution to the LJ energy for this centering atom.
        end
    end

    # Make sure I'm always returning energy per atom.
    if lowercase(LJ.fitTo) == "total"
        return ( (totalEnergy + LJ.offset) * LJ.stdEnergy + LJ.meanEnergy) / crystal.nAtoms
    else
        return (totalEnergy + LJ.offset) * LJ.stdEnergy + LJ.meanEnergy
    end
end



function logNormal(data::DataSet,LJ::LJ,σ::Float64)::Float64
    thesum = 0.0
    for i = 1:data.nData
        thesum += (data.crystals[i].fitEnergy -totalEnergy(data.crystals[i],LJ))^2
    end
    thesum *= - 1/(2 * σ^2)
    
#    @time "thesum" thesum =  -data.nData/2 *log(σ^2) - 1/(2 * σ^2) * sum([(i.energyPerAtomFP - totalEnergy(i,LJ))^2   for i in data.crystals])
    return -data.nData/2 *log(σ^2) + thesum

end

function initializeLJ(settings::Dict)

    order = settings["order"]::Int64
    cutoff = settings["cutoff"]::Float64
    nInteractions = Int(order*(order + 1)/2)
    fitTo = settings["fitTo"]::String
    #params = ones(order,order,2)
    σ = zeros(order,order)
    ϵ = zeros(order,order)
    
    return LJ(order,cutoff,σ,ϵ,1.0,0.0,0.0,fitTo)
end

function initializeLJ(path::String)
    input = YAML.load_file(path,dicttype = Dict{String,Any})

    modelType = string(input["model"]["type"])
    if lowercase(modelType) != "lj"
        error("Model specified is not LJ")
    end
    
    # Get the dataset
    species = String[x for x in input["dataset"]["species"]]
    #species = [x for x in speciesList]
    dataFile = input["dataset"]["file"]
    offset = Float64(input["dataset"]["offset"])
    dset = MatSim.readStructuresIn(dataFile,species,overwriteLatPar = false,offset = offset)
    standardize = Bool(input["dataset"]["standardize"])
    fitTo = String(input["dataset"]["fitTo"])
    #if standardize
    #    MatSim.standardizeData!(dset,offset)
    #end
    # Check to make sure that the specified order matches the dataset
    order = Int(input["model"]["order"])
    if dset.crystals[1].order != order
        error("Order of model not consistent with order of crystals in data set.")
    end

    #offset = Float64(input["dataset"]["offset"])
    # Initialize the LJ model     
    modelDict = input["model"]::Dict{String,Any}
    modelDict["fitTo"] = fitTo
    LJ_model = MatSim.initializeLJ(modelDict)#input["model"])
    # Pre-calculate the distances needed for LJ
    for crystal in dset.crystals
#        MatSim.totalDistances!(crystal,LJ_model)
        MatSim.totalDistances!(crystal,LJ_model)
    end
    # Split data set into training and holdout sets
    nTraining = Int(input["dataset"]["nTraining"])
    #println(nTraining)
    if nTraining > length(dset.crystals)
        error("Data set not big enough for $nTraining training data points")
    end
    trainingSet, holdoutSet = MatSim.getTraining_Holdout_Sets(dset,nTraining,fitTo,standardize)

    # Get everything needed to run Metropolis Hastings.
    metropDict = input["metrop"]::Dict{String,Any}
    LJMetrop = initializeMetrop(input["metrop"],LJ_model) 
    return LJMetrop,trainingSet,holdoutSet,LJ_model,dset
end

function logPost(data,model,metrop,σ)::Float64
    # return 3.0
#    thesum = 0.0
    thesum = MatSim.logNormal(data,model,σ)  + logpdf(metrop.std_Prior,σ)
    for i= 1:model.order, j= i:model.order
        thesum += logpdf(metrop.σ_Priors[i,j],model.σ[i,j])
        thesum += logpdf(metrop.ϵ_Priors[i,j],model.ϵ[i,j])
    end
    return thesum
    #    logPost(data,model,σ) = MatSim.logNormal(data,model,σ) + sum(logpdf.(σ_Priors,model.σ)) + sum(logpdf.(ϵ_Priors,model.ϵ)) + logpdf(std_Prior,σ)  

end

function proposal(μ,σ)
    return Gamma(μ^2/σ^2,σ^2/μ)#proposalDict[lowercase(metrop["proposal"])]
end


function initializeMetrop(metrop::Dict,LJ::LJ)
    order = LJ.order
    nInteractionTypes = Int(order * (order + 1)/2)

    # Check to make sure I have all the priors needed.
    priors = metrop["Priors"]
    if length(priors) - 1 != nInteractionTypes
        error("Number of priors specified is not consistent with the declared order")
    end


    indexDict = Dict("aa"=>[1,1], "ab"=>[1,2], "bb"=>[2,2])
    distDict = Dict("gamma"=> Gamma,"uniform"=>Uniform)
    
    # Get all of the candidate sigmas
    sigs = metrop["candidateSigmas"]
    candSig = zeros(order,order,2)
    ϵ_candSigs = zeros(order,order)
    σ_candSigs = zeros(order,order)
    for i in keys(sigs)
        if lowercase(i) != "sigma"
            ϵ_candSigs[indexDict[i]...] = sigs[i]["epsilon"]
            σ_candSigs[indexDict[i]...] = sigs[i]["sigma"]
#            candSig[indexDict[i]...,1]= sigs[i]["epsilon"]
#            candSig[indexDict[i]...,2]= sigs[i]["sigma"]
        end
    end
    std_candSig = sigs["sigma"]

    # Get all of the staring guesses
    start = metrop["starting"]

    for i in keys(sigs)
        if lowercase(i) != "sigma"
            LJ.ϵ[indexDict[i]...]= start[i]["epsilon"]
            LJ.σ[indexDict[i]...]= start[i]["sigma"]
        end
    end
    std_Guess = start["sigma"]

    # Build the array of priors
    priors = metrop["Priors"]
    paramsPriors = Array{Distribution, 3}(undef, order,order, 2)
    σ_Priors = Array{Gamma{Float64},2}(undef,order,order)
    ϵ_Priors = Array{Gamma{Float64},2}(undef,order,order)
    for i in keys(priors)
        if lowercase(i) != "sigma"
            σ_Priors[indexDict[i]...]= distDict[lowercase(priors[i]["epsilon"]["distribution"])](parse.(Float64,split(priors[i]["epsilon"]["parameters"]))...)
            ϵ_Priors[indexDict[i]...]= distDict[lowercase(priors[i]["sigma"]["distribution"])](parse.(Float64,split(priors[i]["sigma"]["parameters"]))...)
        end
    end
    # The parameters are stored in an order x order x 2 array with only the upper triangle of it actually used. The upper right 
    # entry corresponds to a-b interaction and the lower right corresponds to b-a interactions. But we have to put some kind of a 
    # prior on the lower left entry or the posterior doesn't evaluate right. Tried just copying the prior at 1,2 -> 2,1, but since
    # its a Gamma, evaluating at zero produces infinity. So I switched it to a uniform.  This is a temporary hack though that is
    # specific to a binary case only.  Needs generalized to any order.

#    extras = [a for a in CartesianIndices(LJ.σ) if a[2] < a[1]]
 #   println(extras)
    σ_Priors[2,1] = Gamma(10,0.5)
    ϵ_Priors[2,1] = Gamma(10,0.5)
    std_Prior = distDict[lowercase(priors["sigma"]["distribution"])](parse.(Float64,split(priors["sigma"]["parameters"]))...)

    
    # Define the proposal distribution
    #if lowercase(metrop["proposal"]) == "gamma"
    #end

    nDraws = metrop["nDraws"]
    nBurnIn = metrop["nBurnin"]
    nTotal = nBurnIn + nDraws
    ϵ_accept = zeros(order,order)
    σ_accept = zeros(order,order)
    std_accept = 0.0

    LJ_metrop(nTotal,nBurnIn,ϵ_accept,σ_accept,[std_accept],ϵ_Priors,σ_Priors,std_Prior,ϵ_candSigs,σ_candSigs,std_candSig,std_Guess)
end


function sample_σ!(metrop::LJ_metrop,data::DataSet,LJ::LJ,LJ_next::LJ,std)

    for j = 1:LJ.order, k = j:LJ.order
        #println(LJ.σ,metrop.σ_candSigs)
        cand = rand(proposal(LJ.σ[j,k],metrop.σ_candSigs[j,k]))
        if cand < 0.05
            continue
        end 

        LJ_next.σ[j,k] = cand  # Set sigma equal to the candidate draw

        # Evaluate the posterior at the candidate draw.
        numerator = logPost(data,LJ_next,metrop,std) + log(pdf(proposal(cand,metrop.σ_candSigs[j,k]),LJ.σ[j,k]))
        LJ_next.σ[j,k] = LJ.σ[j,k]  # Set sigma back to the previous draw.
    
        # Evaluate the posterior again.
        denominator =logPost(data,LJ_next,metrop,std)  + log(pdf(proposal(LJ.σ[j,k],metrop.σ_candSigs[j,k]),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0.0 || ((r < 0.0) & (unif < r))  # Accept?
            LJ_next.σ[j,k] = cand
    #        metrop.params_draws[i,j,k,l] = cand   # Yes!
            metrop.σ_accept[j,k] += 1/metrop.nDraws
        end

    end
end



function sample_ϵ!(metrop::LJ_metrop,data::DataSet,LJ::LJ,LJ_next::LJ,std)

    for j = 1:LJ.order, k = j:LJ.order
        cand = rand(proposal(LJ.ϵ[j,k],metrop.ϵ_candSigs[j,k]))
        if cand < 0.05
            continue
        end 

        LJ_next.ϵ[j,k] = cand  # Set epsilon equal to the candidate draw

        # Evaluate the posterior at the candidate draw.
        numerator = logPost(data,LJ_next,metrop,std) + log(pdf(proposal(cand,metrop.ϵ_candSigs[j,k]),LJ.ϵ[j,k]))
        LJ_next.ϵ[j,k] = LJ.ϵ[j,k]  # Set epsilon back to the previous draw.
        
        # Evaluate the posterior again.
        denominator =logPost(data,LJ_next,metrop,std)  + log(pdf(proposal(LJ.ϵ[j,k],metrop.ϵ_candSigs[j,k]),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0 || ((r < 0) & (unif < r))  # Accept?
            LJ_next.ϵ[j,k] = cand
    #        metrop.params_draws[i,j,k,l] = cand   # Yes!
            metrop.ϵ_accept[j,k] += 1/metrop.nDraws
        end

    end
end

function sample_std!(metrop::LJ_metrop,data::DataSet,LJ::LJ,prev_std)

    cand = rand(proposal(prev_std,metrop.std_candSig))
    if cand < 0.05
        return prev_std
    end 
    numerator = logPost(data,LJ,metrop,cand) + log(pdf(proposal(cand,metrop.std_candSig),prev_std))
    denominator =logPost(data,LJ,metrop,prev_std)  + log(pdf(proposal(prev_std,metrop.std_candSig),cand))
    r = numerator - denominator
    unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
    if r >= 0 || ((r < 0) & (unif < r))  # Accept?
        metrop.std_accept[1] += 1/metrop.nDraws
        return cand
        #metrop.σ_draws[i] = cand   # Yes!
    end
    return prev_std

end


function getSamples(metrop::LJ_metrop,data::DataSet,LJ::LJ)
    LJ_next = deepcopy(LJ)
    intDict = Dict(1=>"a",2=>"b")
    #Write the header to the output file
    cDir = pwd()
    println("Opening file  ", joinpath(cDir,"draws.out"))
    system = "System: " * data.title * "\n"
    filename = "draws-LJ." * lstrip(data.title)
    fitTo = "fitTo: " * data.fitTo * "\n"
    io = open(joinpath(cDir,filename),"w")
    standardized = "Standardized: " * string(data.standardized) * "\n"
    mn = data.standardized ? "μ-energy: " * string(data.meanEnergy) * "\n" : "μ-energy: " * string(0.0) * "\n" 
    sn = data.standardized ? "σ-energy: " * string(data.stdEnergy) * "\n" : "σ-energy: " * string(1.0) * "\n"
    offset = data.standardized ? "offset-energy: " * string(data.offset) * "\n" : "offset-energy: " * string(0.0) * "\n"
    cutoff = "cutoff-radius: " * string(LJ.cutoff) * "\n"
    write(io,system)
    write(io,fitTo)
    write(io,standardized)
    write(io,mn)
    write(io,sn)
    write(io,offset)
    write(io,cutoff)
    acceptPosition = mark(io)  # Mark the current position
    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    nSpaces = "%" * string(15 * (2 * nInteractionTypes + 1)- 2) * "s"  # We use 15 spaces per number so let's allocate exactly the right number of spaces for the first line.
    fstring = Printf.Format(nSpaces)
    header = Printf.format(fstring, "\n")
#    header = ""
    for i =1:LJ.order, j = i:LJ.order
        header *= "         ϵ_" *intDict[i] * intDict[j] * "  "
    end

    for i =1:LJ.order, j = i:LJ.order
        header *= "         σ_" *intDict[i] * intDict[j] * "  "
    end
    header *= "         std\n"
#    header = @sprintf "ϵ_aa  ϵ_ab ϵ_bb σ_aa σ_ab σ_bb std\n"
    println(header)
    write(io,header)
    #close(io)
    std_draw = metrop.std
    for i = 1:metrop.nDraws
        sample_σ!(metrop,data,LJ,LJ_next,std_draw)
        sample_ϵ!(metrop,data,LJ,LJ_next,std_draw)
        std_draw = sample_std!(metrop,data,LJ_next,std_draw)
        LJ.σ[:,:] .= LJ_next.σ[:,:]
        LJ.ϵ[:,:] .= LJ_next.ϵ[:,:]
        if i > metrop.nBurnIn
            writeDraw(LJ_next,std_draw,io)
        end
    end

    seek(io,acceptPosition)

    for i= 1:LJ.order, j=i:LJ.order
        printString = @sprintf "%14.2f " metrop.ϵ_accept[i,j] * 100
        write(io,printString)
    end
    for i= 1:LJ.order, j=i:LJ.order
        printString = @sprintf "%14.2f " metrop.σ_accept[i,j] * 100
        write(io,printString)
    end

    printString = @sprintf "%14.2f\n" metrop.std_accept[1] * 100
    write(io,printString)

    println("σ accept")
    display(metrop.σ_accept)
    println("ϵ accept")
    display(metrop.ϵ_accept)
    println("std accept")
    display(metrop.std_accept)
    close(io)

end

function writeDraw(LJ::LJ,std_draw,file)
    
    #io = open(file,"a")
    for i= 1:LJ.order, j=i:LJ.order
        printString = @sprintf "%15.5f" LJ.ϵ[i,j]
        write(file,printString)
    end
    for i= 1:LJ.order, j=i:LJ.order
        printString = @sprintf "%15.5f" LJ.σ[i,j]
        write(file,printString)
    end
    
    printString = @sprintf "%15.5f\n" std_draw
    write(file,printString)
end

