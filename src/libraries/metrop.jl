#module metrop

# Initialize directly from the YAML file.
function initializeModel(path::String)
    input = YAML.load_file(path)
    
    species = input["dataset"]["species"]
    order = input["model"]["order"]
#    nParams = Int(order * (order + 1)/2) * 2
#    priors = input["metrop"]["Priors"]
#    if length(priors) - 1 != nParams
#        error("Number of priors specified is not consistent with the declared order")
#    end
    dataFolder = dirname(input["dataset"]["file"])
    dataFile = basename(input["dataset"]["file"])
    dset = MatSim.readStructuresIn(dataFolder,dataFile,species,overwriteLatPar = false)
    if dset.crystals[1].order != order
        println(order)
        error("Order of model not consistent with order of crystals in data set.")
    end
    offset = input["dataset"]["offset"]
    #display(dset)
    ## Standardize the data..
    #if input["dataset"]["standardize"]
    #    MatSim.standardizeData!(dset,offset)
    #end
    # Model dictionary
    modelDict = Dict("LJ"=> MatSim.initializeLJ)
    energyModel = modelDict[input["model"]["type"]](input["model"])
    ## Pre-calculate the distances needed for LJ
#    cutoff = input["model"]["cutoff"]
    for crystal in dset.crystals
        crystal.ljvals .= MatSim.totalDistances(crystal,energyModel)
    end
    # Split data set into training and holdout sets
    nTraining = input["dataset"]["nTraining"]
    if nTraining > length(dset.crystals)
        error("Data set not big enough for $nTraining training data points")
    end
    trainingSet, holdoutSet = MatSim.getTraining_Holdout_Sets(dset,nTraining,standardize = input["dataset"]["standardize"],offset = offset)

    thisMetrop = initializeMetrop(input["metrop"],energyModel) 
    return thisMetrop,trainingSet,holdoutSet,energyModel#model(energyModel,thisMetrop,trainingSet,holdoutSet)
end
# Reads the settings dictionary and builds all the settings variables.
# Then initializes a LJ object for use in a MH algorithm
function initializeMetrop(metropSettings::Dict,LJ::LJ)
    order = LJ.order
    nInteractionTypes = Int(order * (order + 1)/2)

    priors = metropSettings["Priors"]
    if length(priors) - 1 != nInteractionTypes
        error("Number of priors specified is not consistent with the declared order")
    end
    distDict = Dict("gamma"=> Gamma)
    sigs = metropSettings["candidateSigmas"]
    nDraws = metropSettings["nDraws"]
    nBurnIn = metropSettings["nBurnin"]
    candSig = [sigs["aa"]["sigma"] sigs["aa"]["epsilon"]
               sigs["ab"]["sigma"] sigs["ab"]["epsilon"]
               sigs["bb"]["sigma"] sigs["bb"]["epsilon"]
              ]
    candSig_sig = sigs["sigma"]

    start = metropSettings["starting"]
    display(start)
    μGuess = [start["aa"]["sigma"] start["aa"]["epsilon"]
              start["ab"]["sigma"] start["ab"]["epsilon"]
              start["bb"]["sigma"] start["bb"]["epsilon"]
             ]

    σGuess = start["sigma"]


    priors = metropSettings["Priors"]
    #order = Int((length(priors)-1)/2)
    μprior = Array{Distribution, 2}(undef, nInteractionTypes, 2)
    μprior[1,1] =  distDict[lowercase(priors["aa"]["sigma"]["distribution"])](parse.(Int,split(priors["aa"]["sigma"]["parameters"]))...)
    μprior[1,2] =  distDict[lowercase(priors["aa"]["epsilon"]["distribution"])](parse.(Int,split(priors["aa"]["epsilon"]["parameters"]))...)
    μprior[2,1] =  distDict[lowercase(priors["ab"]["sigma"]["distribution"])](parse.(Int,split(priors["ab"]["sigma"]["parameters"]))...)
    μprior[2,2] =  distDict[lowercase(priors["ab"]["epsilon"]["distribution"])](parse.(Int,split(priors["ab"]["epsilon"]["parameters"]))...)
    μprior[3,1] =  distDict[lowercase(priors["bb"]["sigma"]["distribution"])](parse.(Int,split(priors["bb"]["sigma"]["parameters"]))...)
    μprior[3,2] =  distDict[lowercase(priors["bb"]["epsilon"]["distribution"])](parse.(Int,split(priors["bb"]["epsilon"]["parameters"]))...)
    σprior = distDict[lowercase(priors["sigma"]["distribution"])](parse.(Int,split(priors["sigma"]["parameters"]))...)
    logPost(data,energyModel,σ) = MatSim.loglik(data,energyModel,σ) +  sum(logpdf.(μprior,energyModel.params)) + logpdf(σprior,σ)
    
    if lowercase(metropSettings["proposal"]) == "gamma"
        proposal(μ,σ) = Gamma(μ^2/σ^2,σ^2/μ)#proposalDict[lowercase(metropSettings["proposal"])]
    end



    nTotal = nBurnIn + nDraws
    nParams = size(candSig)
    μ_draws = zeros(nTotal,nParams...)
    μ_draws[1,:,:] .= μGuess
    σ_draws = zeros(nTotal)
    σ_draws[1] = σGuess
    μ_accept = zeros(nParams...)
    σ_accept = Float64[0]

    metrop(nTotal,nBurnIn,μ_draws,σ_draws,candSig,candSig_sig,μGuess,σGuess,μ_accept,σ_accept, proposal,logPost)
end



function getSamples(metropSettings::metrop,data::DataSet,energyModel::LJ)
    nParams = size(metropSettings.candSig_μ)
#    metrop.model.params .= zeros(nParams...)
    
#    drawsWithCand = zeros(nParams...)
    for i = 2:metropSettings.nDraws
        metropSettings.μ_draws[i,:,:] .= metropSettings.μ_draws[i-1,:,:]  # Set the next draw to be equal to the previous.  I
        energyModel.params .= metropSettings.μ_draws[i,:,:]  # Need to assemble the vector of parameters with the candidate draw inserted at the right place.
        metropSettings.σ_draws[i] = metropSettings.σ_draws[i-1]
        for j = 1: nParams[1], k = 1:nParams[2]
            cand = rand(metropSettings.proposal(metropSettings.μ_draws[i,j,k],metropSettings.candSig_μ[j,k]))
            if cand < 0.05
                continue
            end 

            energyModel.params[j,k] = cand
            numerator = metropSettings.logpost(data,energyModel,metropSettings.σ_draws[i]) + log(pdf(metropSettings.proposal(cand,metropSettings.candSig_μ[j,k]),metropSettings.μ_draws[i,j,k]))
            energyModel.params[j,k] = metropSettings.μ_draws[i,j,k]
            
            denominator =metropSettings.logpost(data,energyModel,metropSettings.σ_draws[i])  + log(pdf(metropSettings.proposal(metropSettings.μ_draws[i,j,k],metropSettings.candSig_μ[j,k]),cand))

            r = numerator - denominator
            unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
            if r >= 0 || ((r < 0) & (unif < r))  # Accept?
                energyModel.params[j,k] = cand
                metropSettings.μ_draws[i,j,k] = cand   # Yes!
                metropSettings.μ_accept[j,k] += 1/metropSettings.nDraws
            end
        end

        energyModel.params .= metropSettings.μ_draws[i,:,:]  # Get the updated mu values to use when getting draws for sigma

        ## Now get sigma draws...

        cand = rand(metropSettings.proposal(metropSettings.σ_draws[i],metropSettings.candSig_σ))
        if cand < 0.05
            continue
        end 
        numerator = metropSettings.logpost(data,energyModel,cand) + log(pdf(metropSettings.proposal(cand,metropSettings.candSig_σ),metropSettings.σ_draws[i]))
        denominator =metropSettings.logpost(data,energyModel,metropSettings.σ_draws[i])  + log(pdf(metropSettings.proposal(metropSettings.σ_draws[i],metropSettings.candSig_σ),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0 || ((r < 0) & (unif < r))  # Accept?
            metropSettings.σ_draws[i] = cand   # Yes!
            metropSettings.σ_accept[1] += 1/metropSettings.nDraws
        end

    end     
    return metropSettings
end

#end



#function getSamples(metrop::Metrop)
#    drawsWithCand = zeros(metrop.nParams)
#    for i = 1:metrop.nDraws
#        metrop.draws[i,:] .= metrop.draws[i-1,:]  # Set the next draw to be equal to the previous.  I
#        for j = 1: metrop.nParams
#            cand = rand(Normal(draws[i,j],metrop.candSig[j]))  # Get a candidate draw.  Draw from distribution with previous draw
#            drawsWithCand[:] .= metrop.draws[:]  # Need to assemble the vector of parameters with the candidate draw inserted at the right place.
#            drawsWithCand[j] = cand
#            r = metrop.logpost(data,drawsWithCand...) - logga(data,metrop.draws[i,:]...)  #Ratio between last draw with 
#            unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
#            if r >= 0 || (r < 0 && unif < r)  # Accept?
#                metrop.draws[i,j] .= cand   # Yes!
#                metrop.accept[j] .+= 1/metrop.nDraws
#            end
#        end
#    end     
#    println("acceptance rates", metrop.accept * 100)
#    return metrop
#end


#end

#using Pkg
#Pkg.add("Plots")
#using Distributions
#using DelimitedFiles
#using Plots
#using SpecialFunctions
#
#
#cd("/Users/legoses/Downloads/")
#data = readdlm("projectile.csv",',',Float64) # Read in the data
#μ = mean(data[:,1])
#σ = std(data[:,1])
#data[:,1] .= (data[:,1] .- mean(data[:,1]))./std(data[:,1])
#data[:,2] .= (data[:,2] .- mean(data[:,2]))./std(data[:,2])
#draws = bigLoop(data)
#
#
#a = (draws[:,1] .- μ)/σ
#histogram(a)
#plot(draws[:,3])
#plot(data[:,1],data[:,2],seriestype = :scatter)
#x = 0:0.01:10
#x = (x .- mean(x))/std(x)
#y = [draws[i,1] .* x.^2 .+ draws[i,2] .* x .+ draws[i,3] for i in 20000:1000:nDraws]
#display(draws[:,1])
#plot(x,y,legend=false)
#histogram(draws[10000:end,4],bins = 200)
#histogram2d(draws[10000:end,1], draws[10000:end,2])
#
#maximum(draws[2000:,2])
#minimum(draws[:,2])
#xPred = 175
#yPred = [rand(Normal(draws[i,1] * xPred^2 + draws[i,2]*xPred + draws[i,3],sqrt(21))) for i in 2000:100000]
#histogram(yPred)
#
#
#a = 0:0.01:3
#y = [logga(data,x,4,15) for x in a]
#maximum(y)
#plot(a,exp.(y))
#
#
#