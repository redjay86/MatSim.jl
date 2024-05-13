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
                indices = sort([iNeighbor,centerType])
                if all([i,j,k] .== 0 ) 
                    totalEnergy -=  4 * LJ.params[indices...,1] * 1/2 * LJ.params[indices...,2]^6/r^6
                    totalEnergy +=  4 * LJ.params[indices...,1] * 1/2 * LJ.params[indices...,2]^12/r^12
                else 
                    totalEnergy -= 4 * LJ.params[indices...,1] * LJ.params[indices...,2]^6/r^6
                    totalEnergy += 4 * LJ.params[indices...,1] * LJ.params[indices...,2]^12/r^12
                end
            end
        end
    end
    return totalEnergy
end

function singleAtomDistances!(distMat::Array{Float64,3},crystal::Crystal,LJ::LJ,centerAtom::Vector{Float64}, centerType:: Integer, loopBounds::Vector{Int64})
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
                # The LJ parameters are stored in the upper triangular portion of a matrix
                # The bottom triangle is redundant.. interaction between a and b is equivalent
                # to interaction between b and a.  So I sort the indices here so that the bottom
                # triangle of the matrix never gets updated, only the upper right.
                indices = sort([iNeighbor,centerType])
                if all([i,j,k] .== 0 ) 
                    distMat[indices...,1] +=  4 * 1.0/2.0 * 1.0/r^6
                    distMat[indices...,2] +=  4* 1.0/2.0 * 1.0/r^12
                else 
                    distMat[indices...,1] += 4 * 1.0/r^6
                    distMat[indices...,2] += 4 * 1.0/r^12
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
        totalEnergy = sum(-LJ.params[:,:,1] .* LJ.params[:,:,2] .^6 .* crystal.ljvals[:,:,1] + LJ.params[:,:,1] .* LJ.params[:,:,2] .^12 .* crystal.ljvals[:,:,2] ) 
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



function logNormal(data::DataSet,LJ::LJ,σ::Float64)
    thesum =  -data.nData/2 *log(σ^2) - 1/(2 * σ^2) * sum([(i.energyFP - totalEnergy(i,LJ))^2   for i in data.crystals])
    return thesum

end

function initializeLJ(settings::Dict)

    order = settings["order"]
    cutoff = settings["cutoff"]
    nInteractions = Int(order*(order + 1)/2)
    params = ones(order,order,2)
    return LJ(order,params,cutoff)
end

function initializeLJ(path::String)
    input = YAML.load_file(path)
    if lowercase(input["model"]["type"]) != "lj"
        error("Model specified is not LJ")
    end
    
    # Get the dataset
    species = input["dataset"]["species"]
    dataFolder = dirname(input["dataset"]["file"])
    dataFile = basename(input["dataset"]["file"])
    dset = MatSim.readStructuresIn(dataFolder,dataFile,species,overwriteLatPar = false)
    if input["dataset"]["standardize"]
        MatSim.standardizeData!(dset,input["dataset"]["offset"])
    end
    # Check to make sure that the specified order matches the dataset
    order = input["model"]["order"]
    if dset.crystals[1].order != order
        println(order)
        error("Order of model not consistent with order of crystals in data set.")
    end

    offset = input["dataset"]["offset"]
    # Initialize the LJ model     
    LJ = MatSim.initializeLJ(input["model"])
    ## Pre-calculate the distances needed for LJ
    for crystal in dset.crystals
        crystal.ljvals .= MatSim.totalDistances(crystal,LJ)
    end
    # Split data set into training and holdout sets
    nTraining = input["dataset"]["nTraining"]
    if nTraining > length(dset.crystals)
        error("Data set not big enough for $nTraining training data points")
    end
    trainingSet, holdoutSet = MatSim.getTraining_Holdout_Sets(dset,nTraining)

    # Get everything needed to run Metropolis Hastings.
    LJMetrop = initializeMetrop(input["metrop"],LJ) 
    return LJMetrop,trainingSet,holdoutSet,LJ,dset
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
    for i in keys(sigs)
        if lowercase(i) != "sigma"
            candSig[indexDict[i]...,1]= sigs[i]["epsilon"]
            candSig[indexDict[i]...,2]= sigs[i]["sigma"]
        end
    end
    candSig_sig = sigs["sigma"]

    # Get all of the staring guesses
    start = metrop["starting"]
    params_Guess = zeros(order,order,2)
    for i in keys(sigs)
        if lowercase(i) != "sigma"
            params_Guess[indexDict[i]...,1]= start[i]["epsilon"]
            params_Guess[indexDict[i]...,2]= start[i]["sigma"]
        end
    end
    σGuess = start["sigma"]


    # Build the array of priors
    priors = metrop["Priors"]
    paramsPriors = Array{Distribution, 3}(undef, order,order, 2)
    for i in keys(priors)
        println(i)
        if lowercase(i) != "sigma"
            paramsPriors[indexDict[i]...,1]= distDict[lowercase(priors[i]["epsilon"]["distribution"])](parse.(Float64,split(priors[i]["epsilon"]["parameters"]))...)
            paramsPriors[indexDict[i]...,2]= distDict[lowercase(priors[i]["sigma"]["distribution"])](parse.(Float64,split(priors[i]["sigma"]["parameters"]))...)
        end
    end
    # The parameters are stored in an order x order x 2 array with only the upper triangle of it actually used. The upper right 
    # entry corresponds to a-b interaction and the lower right corresponds to b-a interactions. But we have to put some kind of a 
    # prior on the lower left entry or the posterior doesn't evaluate right. Tried just copying the prior at 1,2 -> 2,1, but since
    # its a Gamma, evaluating at zero produces infinity. So I switched it to a uniform.  This is a temporary hack though that is
    # specific to a binary case only.  Needs generalized to any order.

    extras = [a for a in CartesianIndices(LJ.params) if a[2] < a[1]]
    paramsPriors[extras] .= Uniform(-5,5)
    σprior = distDict[lowercase(priors["sigma"]["distribution"])](parse.(Float64,split(priors["sigma"]["parameters"]))...)

    # Construct the posterior
    logPost(data,model,σ) = MatSim.logNormal(data,model,σ) + sum(logpdf.(paramsPriors,model.params)) + logpdf(σprior,σ)  
    
    # Define the proposal distribution
    if lowercase(metrop["proposal"]) == "gamma"
        proposal(μ,σ) = Gamma(μ^2/σ^2,σ^2/μ)#proposalDict[lowercase(metrop["proposal"])]
    end

    nDraws = metrop["nDraws"]
    nBurnIn = metrop["nBurnin"]
    nTotal = nBurnIn + nDraws
#    nParams = size(candSig)
    params_draws = zeros(nTotal,order,order,2)
    params_draws[1,:,:,:] .= params_Guess
    σ_draws = zeros(nTotal)
    σ_draws[1] = σGuess
    params_accept = zeros(order,order,2)
    σ_accept = Float64[0]

    LJ_metrop(nTotal,nBurnIn,params_draws,σ_draws,candSig,candSig_sig,params_Guess,σGuess,params_accept,σ_accept, proposal,logPost,paramsPriors,σprior)
end


function displayResults(results::LJ_metrop,LJ::LJ,trainingSet::DataSet,holdoutSet::DataSet,meanEnergy::Float64, stdEnergy::Float64,offset::Float64)
    #histogram(results.σ_draws)
    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)
    if size(LJ.params)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end
    # We don't use the lower left triangle of the matrix of parameters, so let's get the indices right now and sort them so the plots
    # are arranged correctly.
    keeps = [a for a in CartesianIndices(LJ.params) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    x = 0:0.01:3
    hists = [histogram(results.params_draws[results.nBurnIn:end,a],bins = 100,normalize = :pdf,annotations = ((0.5,0.5),@sprintf("Acceptance Rate: %5.1f %%",results.params_accept[a]*100))) for a in keeps]
    σ_hist = histogram(results.σ_draws,bins = 100,normalize = :pdf,annotations = ((0.5,0.5),@sprintf("Acceptance Rate: %5.1f %%",results.σ_accept[1]*100)))
    plot!(x,pdf(results.σ_Prior,x),lw =6,lc = :red)
    
    combs = [[[i,j,k] for k = 1:2] for i = 1:LJ.order for j = i:LJ.order ]

    hist2ds = [histogram2d(results.params_draws[results.nBurnIn:end,x[1]...],results.params_draws[results.nBurnIn:end,x[2]...]) for x in combs  ]
    r = @layout [ a; b; c]
    hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = 1)
    l = @layout [grid(nInteractionTypes,2)]
    resplot = plot(hists...,layout = l,size = (2250,1250),legend=false)
    plot!(x,[pdf(results.μ_Priors[a],x) for a in keeps],lw=6,lc = :red)
    g = @layout [grid(nInteractionTypes,2)]
    tracePlots = plot([results.params_draws[results.nBurnIn:end,a] for a in keeps],layout = g, size = (2000,1000),legend = false)
    #tracePlots = plot(trace...,layout = l,size = (2000,1000),legend=false) 
    
    # Predict on the holdout set to evaluate quality of model.
#    stdEnergy = trainingSet.stdEnergy
#    meanEnergy = trainingSet.meanEnergy
#    offset = trainingSet.offset
    trueVals = zeros(Float64,length(holdoutSet.crystals))
    predictVals = zeros(Float64,length(holdoutSet.crystals))
    predictUnc = zeros(Float64,length(holdoutSet.crystals))
    rmsError = zeros(Float64,length(holdoutSet.crystals))
    for j = 1:length(holdoutSet.crystals)
        trueVals[j] = (holdoutSet.crystals[j].energyFP + offset) * stdEnergy + meanEnergy 
#        trueVals[j] = holdoutSet.crystals[j].energyFP + offset
        overDraws = zeros(Float64,results.nDraws - results.nBurnIn)
        for i = 1:results.nDraws - results.nBurnIn
            LJ.params .= results.params_draws[i+results.nBurnIn,:,:,:]
#            overDraws[i] = MatSim.totalEnergy(holdoutSet.crystals[j],LJ)
            overDraws[i] = (MatSim.totalEnergy(holdoutSet.crystals[j],LJ) + offset) * stdEnergy + meanEnergy
        end
        predictVals[j] = mean(overDraws)
        predictUnc[j] = std(overDraws)
    end
    rmsError = sqrt(mean( (trueVals .- predictVals).^2 ))
    upper = maximum(trueVals)
    lower = minimum(trueVals)
    x = 1.15*lower:0.05:0.85*upper
    
    # Plot predicted vs true energies. Slope=1 line is a perfect fit.
    myp = plot(predictVals,trueVals,seriestype = :scatter,xerror = predictUnc,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "Predicted Energy (eVs/atom)",legend=false)
    plot!(x,x,lw = 5)
    
 #   using Printf
    tag = @sprintf("RMS Error: %5.2f",rmsError)
    annotate!((0.75,0.25),tag)
    return resplot,myp,σ_hist,tracePlots,hist2dplots
        
end


function getSamples(metrop::LJ_metrop,data::DataSet,LJ::LJ)
    nParams = size(metrop.candSig_params)
#    metrop.model.params .= zeros(nParams...)
    order = LJ.order
#    drawsWithCand = zeros(nParams...)
    for i = 2:metrop.nDraws
       # println(i)
        metrop.params_draws[i,:,:,:] .= metrop.params_draws[i-1,:,:,:]  # Set the next draw to be equal to the previous.  I
        LJ.params .= metrop.params_draws[i,:,:,:]  # Need to assemble the vector of parameters with the candidate draw inserted at the right place.
        metrop.σ_draws[i] = metrop.σ_draws[i-1]
        for j = 1:order, k = j:order, l = 1:2
            cand = rand(metrop.proposal(metrop.params_draws[i,j,k,l],metrop.candSig_params[j,k,l]))
            if cand < 0.05
                continue
            end 

            LJ.params[j,k,l] = cand
            numerator = metrop.logpost(data,LJ,metrop.σ_draws[i]) + log(pdf(metrop.proposal(cand,metrop.candSig_params[j,k,l]),metrop.params_draws[i,j,k,l]))
            LJ.params[j,k,l] = metrop.params_draws[i,j,k,l]
            
            denominator =metrop.logpost(data,LJ,metrop.σ_draws[i])  + log(pdf(metrop.proposal(metrop.params_draws[i,j,k,l],metrop.candSig_params[j,k,l]),cand))
            r = numerator - denominator
            unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
            if r >= 0 || ((r < 0) & (unif < r))  # Accept?
                LJ.params[j,k,l] = cand
                metrop.params_draws[i,j,k,l] = cand   # Yes!
                metrop.params_accept[j,k,l] += 1/metrop.nDraws
            end
        end

        #LJ.params .= metrop.params_draws[i,:,:,:]  # Get the updated mu values to use when getting draws for sigma

        ## Now get sigma draws...

        cand = rand(metrop.proposal(metrop.σ_draws[i],metrop.candSig_σ))
        if cand < 0.05
            continue
        end 
        numerator = metrop.logpost(data,LJ,cand) + log(pdf(metrop.proposal(cand,metrop.candSig_σ),metrop.σ_draws[i]))
        denominator =metrop.logpost(data,LJ,metrop.σ_draws[i])  + log(pdf(metrop.proposal(metrop.σ_draws[i],metrop.candSig_σ),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0 || ((r < 0) & (unif < r))  # Accept?
            metrop.σ_draws[i] = cand   # Yes!
            metrop.σ_accept[1] += 1/metrop.nDraws
        end

    end     
    return metrop
end

#end