

function σ_hists(results::LJ_metrop,LJ::LJ)

    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),skipstart = 2)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    σ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    for i = 1:LJ.order, j = i:LJ.order
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + (i - 1) * LJ.order + j])
    end
    if size(LJ.σ)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end

    keeps = [a for a in CartesianIndices(LJ.σ) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3
    σ_hists = [histogram(σ_draws[:,a],bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],results.σ_accept[a]*100),6))) for a in keeps]

    σ = plot(σ_hists...)
    plot!(x,[pdf(results.σ_Priors[a],x) for a in keeps],lw=6,lc = :red)

    return σ
end


function ϵ_hists(results::LJ_metrop,LJ::LJ)

    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),skipstart = 2)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    ϵ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    for i = 1:LJ.order, j = i:LJ.order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,(i - 1) * LJ.order + j])
    end
    if size(LJ.σ)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end

    keeps = [a for a in CartesianIndices(LJ.σ) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3
    ϵ_hists = [histogram(ϵ_draws[:,a],bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],results.ϵ_accept[a]*100),6))) for a in keeps]

    ϵ = plot(ϵ_hists...)
    plot!(x,[pdf(results.ϵ_Priors[a],x) for a in keeps],lw=6,lc = :red)


    return ϵ
end

function std_hist(results::LJ_metrop,LJ::LJ)

    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),skipstart = 2)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    std_draws = convert.(Float64,data[:,end])
    if size(LJ.σ)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end


    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3


    std_hist = histogram(std_draws,bins = 100,normalize = :pdf,annotations = ((0.5,0.5),@sprintf("Acceptance Rate: %5.1f %%",results.std_accept[1]*100)))
    plot!(x,pdf(results.std_Prior,x),lw =6,lc = :red)
    return std_hist
end


function LJAverages(results,order,cutoff)
    nInteractionTypes = Int(order * (order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),Float64;skipstart = 2)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    ϵ_mean = zeros(order,order)
    σ_mean = zeros(order,order)
    for i = 1:order, j = i:order
        display(ϵ_mean)
        display(σ_mean)
        for k = 1:size(data)[1]  # Loop over all the draws
    #        println(data[k,:])
            ϵ_mean[i,j] += convert.(Float64,data[k,(i - 1) * order + j])/results.nDraws
            σ_mean[i,j] += convert.(Float64,data[k,nInteractionTypes + (i - 1) * order + j])/results.nDraws
        end
    end

    return LJ(order, cutoff,σ_mean,ϵ_mean)
end
function predPlot(results::LJ_metrop,LJ::LJ,trainingSet::DataSet,holdoutSet::DataSet,meanEnergy::Float64, stdEnergy::Float64,offset::Float64)
#histogram(results.σ_draws)
    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),Float64;skipstart = 2)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    ϵ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    σ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    for i = 1:LJ.order, j = i:LJ.order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,(i - 1) * LJ.order + j])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + (i - 1) * LJ.order + j])
    end
    if size(LJ.σ)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end

    trueVals = zeros(Float64,length(holdoutSet.crystals))
    predictVals = zeros(Float64,length(holdoutSet.crystals))
    predictUnc = zeros(Float64,length(holdoutSet.crystals))
    rmsError = zeros(Float64,length(holdoutSet.crystals))
    for j = 1:length(holdoutSet.crystals)
        trueVals[j] = (holdoutSet.crystals[j].energyFP + offset) * stdEnergy + meanEnergy 
        overDraws = zeros(Float64,results.nDraws - results.nBurnIn)
        for i = 1:results.nDraws - results.nBurnIn
            LJ.ϵ[:,:] .= ϵ_draws[i,:,:]
            LJ.σ[:,:] .= σ_draws[i,:,:]
            overDraws[i] = (MatSim.totalEnergy(holdoutSet.crystals[j],LJ) + offset) * stdEnergy + meanEnergy
        end
        predictVals[j] = mean(overDraws)
        predictUnc[j] = std(overDraws)
    end
    percentError = (predictVals .- trueVals)./trueVals
    percentHist = histogram(percentError,bins = 100,normalize = :pdf,title = "percent error")
    rmsError = sqrt(mean( (trueVals .- predictVals).^2 ))


    upper = maximum(trueVals)
    lower = minimum(trueVals)
    x = 1.15*lower:0.05:0.85*upper
    r = @layout [grid(2,1)]

    # Plot predicted vs true energies. Slope=1 line is a perfect fit.
    myp = plot(predictVals,trueVals,seriestype = :scatter,xerror = predictUnc,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "Predicted Energy (eVs/atom)",legend=false)
    tag = @sprintf("RMS Error: %8.4f",rmsError)
    annotate!((0.75,0.25),tag)
    final = plot(myp,percentHist,layout = r)
    plot!(x,x,lw = 5)
    #using Printf
    return final
end

function tracePlots(results,LJ)

    #histogram(results.σ_draws)
    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),skipstart = 2)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    ϵ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    σ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    std_draws = convert.(Float64,data[:,end])
    for i = 1:LJ.order, j = i:LJ.order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,(i - 1) * LJ.order + j])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + (i - 1) * LJ.order + j])
    end
    if size(LJ.σ)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end
    # We don't use the lower left triangle of the matrix of parameters, so let's get the indices right now and sort them so the plots
    # are arranged correctly.
    keeps = [a for a in CartesianIndices(LJ.σ) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    g = @layout [grid(LJ.order,LJ.order)]
    intDict = Dict(1=>"a",2=>"b")
    #,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],results.ϵ_accept[a]*100),6))
    tracePlots = plot([ϵ_draws[results.nBurnIn:end,a] for a in keeps]
,layout = g, size = (2000,1000),legend = false)
    plot!( [σ_draws[results.nBurnIn:end,a] for a in keeps],layout = g, size = (2000,1000),legend = false)
    return tracePlots
end



function hists2d(results::LJ_metrop,LJ::LJ,type; ar = 1.0)
    #histogram(results.σ_draws)
    nInteractionTypes = Int(LJ.order * (LJ.order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    #Read the draws from file
    data = readdlm(joinpath(cDir,"draws.out"),skipstart = 2,Float64)
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    ϵ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    σ_draws = zeros(results.nDraws-results.nBurnIn,LJ.order,LJ.order)
    std_draws = convert.(Float64,data[:,end])
    for i = 1:LJ.order, j = i:LJ.order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,(i - 1) * LJ.order + j])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + (i - 1) * LJ.order + j])
    end
    if size(LJ.σ)[1] != LJ.order
        error("Number of interaction types not matching up with specified order")
    end

    intDict = Dict(1=>"a",2=>"b")
    if type == "σ-σ" || type == "ϵ-ϵ"
        combs = collect(multiset_combinations(1:nInteractionTypes,2))
        elem = [[i,j] for i = 1:LJ.order for j = i:LJ.order]
        final = [[elem[i[1]],elem[i[2]]] for i in combs]

        if type == "σ-σ"
            hist2ds = [histogram2d(σ_draws[:,x[1]...],σ_draws[:,x[2]...],xlabel = @sprintf("σ-%s%s ",intDict[x[1][1]],intDict[x[1][2]]),ylabel = @sprintf("σ-%s%s",intDict[x[2][1]],intDict[x[2][2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(σ_draws[:,x[1]...]),1.1* maximum(σ_draws[:,x[1]...])),ylim = (0.9*minimum(σ_draws[:,x[2]...]),1.1* maximum(σ_draws[:,x[2]...]))) for x in final  ]
            r = @layout [grid(length(combs),1)] 
            hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar, size = (2000,1000))
        else
            hist2ds = [histogram2d(ϵ_draws[:,x[1]...],ϵ_draws[:,x[2]...],xlabel = @sprintf("ϵ-%s%s ",intDict[x[1][1]],intDict[x[1][2]]),ylabel = @sprintf("ϵ-%s%s",intDict[x[2][1]],intDict[x[2][2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(ϵ_draws[:,x[1]...]),1.1* maximum(ϵ_draws[:,x[1]...])),ylim = (0.9*minimum(ϵ_draws[:,x[2]...]),1.1* maximum(ϵ_draws[:,x[2]...]))) for x in final  ]
            r = @layout [grid(length(combs),1)] 
            hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar, size = (2000,1000))
        end
    else
        elem = [[i,j] for i = 1:LJ.order for j = i:LJ.order]
        hist2ds = [histogram2d(σ_draws[:,y...],ϵ_draws[:,x...],xlabel = @sprintf("σ-%s%s ",intDict[x[1]],intDict[x[2]]),ylabel = @sprintf("ϵ-%s%s",intDict[y[1]],intDict[y[2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(σ_draws[:,y...]),1.1* maximum(σ_draws[:,y...])),ylim = (0.9*minimum(ϵ_draws[:,x...]),1.1* maximum(ϵ_draws[:,x...]))) for x in elem for y in elem ]
        r = @layout [grid(nInteractionTypes,nInteractionTypes)] 
        hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar,colorbar=false, size = (2000,1000))

    end
    return hist2dplots
        
end
