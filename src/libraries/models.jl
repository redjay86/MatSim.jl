# Initialize directly from the YAML file.
function initializeModel(path::String)
    input = YAML.load_file(path)
    
    species = input["dataset"]["species"]
    dataFolder = dirname(input["dataset"]["file"])
    dataFile = basename(input["dataset"]["file"])
    dset = MatSim.readStructuresIn(dataFolder,dataFile,species,overwriteLatPar = false)

    order = input["model"]["order"]
    if dset.crystals[1].order != order
        println(order)
        error("Order of model not consistent with order of crystals in data set.")
    end
    offset = input["dataset"]["offset"]
    # Model dictionary
    modelDict = Dict("LJ"=> MatSim.initializeLJ)
    energyModel = modelDict[input["model"]["type"]](input["model"])
    ## Pre-calculate the distances needed for LJ
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
    return thisMetrop,trainingSet,holdoutSet,energyModel
end

function displayResults(results::metrop,model::LJ,trainingSet::DataSet,holdoutSet::DataSet;offset = 0.0)
    #histogram(results.σ_draws)
    nInteractionTypes = Int(model.order * (model.order + 1)/2)
    if size(model.params)[1] != nInteractionTypes
        error("Number of interaction types not matching up with specified order")
    end
    locations = Array{Vector{Float64}}(undef, nInteractionTypes,2)
    drawMean = mean.(eachslice(results.μ_draws,dims = (2,3)))
    for a in CartesianIndices(drawMean)
        locations[a] = [drawMean[a],1.5]
    end

    display(results.μ_accept)
    hists = [histogram(results.μ_draws[results.nBurnIn:end,a],bins = 100,normalize = :pdf,annotations = (locations[a,:][1]...,@sprintf("Acceptance Rate: %5.1f %%",results.μ_accept[a]*100))) for a in CartesianIndices(model.params)]
    
    tracePlots = [plot(results.μ_draws[results.nBurnIn:end,a],serietype= :scatter)  for a in CartesianIndices(model.params)]
    l = @layout [grid(nInteractionTypes,2)]#[a c; b d; e f]
#    tag = [@sprintf("Acceptance Rate: %5.1f %%",results.μ_accept[a]*100) for a in CartesianIndices(model.params)]
    resplot = plot(hists[[1,4,2,5,3,6]]...,layout = l,size = (2250,1250),legend=false)

    
    # Predict on the holdout set to evaluate quality of model.
    stdEnergy = trainingSet.stdEnergy
    meanEnergy = trainingSet.meanEnergy
    trueVals = zeros(Float64,length(holdoutSet.crystals))
    predictVals = zeros(Float64,length(holdoutSet.crystals))
    predictUnc = zeros(Float64,length(holdoutSet.crystals))
    rmsError = zeros(Float64,length(holdoutSet.crystals))
    for j = 1:length(holdoutSet.crystals)
        trueVals[j] = (holdoutSet.crystals[j].energyFP + offset) * stdEnergy + meanEnergy 
        overDraws = zeros(Float64,results.nDraws - results.nBurnIn)
        for i = 1:results.nDraws - results.nBurnIn
            model.params .= results.μ_draws[i,:,:]
            overDraws[i] = (MatSim.totalEnergy(holdoutSet.crystals[j],model) + offset) *stdEnergy + meanEnergy
        end
        predictVals[j] = mean(overDraws)
        predictUnc[j] = std(overDraws)
    end
    rmsError = sqrt(mean( (trueVals .- predictVals).^2 ))
    x = -7:0.05:-3
    
    # Plot predicted vs true energies. Slope=1 line is a perfect fit.
    myp = plot(predictVals,trueVals,seriestype = :scatter,xerror = predictUnc,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "Predicted Energy (eVs/atom)",legend=false)
    plot!(x,x,lw = 5)
    
 #   using Printf
    tag = @sprintf("RMS Error: %5.2f",rmsError)
    annotate!(-6,-3.2,tag)
    return resplot,myp
        
end