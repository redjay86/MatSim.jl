#module matTypes

#using LinearAlgebra
#using StaticArrays


abstract type energyModel end
# Crystal Structure
mutable struct Crystal
    title::String
    latpar::Float64
    lVecs:: Matrix{Float64}
    nType:: Vector{Int64} #Number of each type of atom 
    aType:: Vector{Int64} # Integers representing the type for each basis atom
    nAtoms::Int64  # Total number of atoms in the unit cell
    coordSys::Vector{String} # 'D' for Direct or 'C' for cartesian
    atomicBasis:: Vector{Vector{Vector{Float64}}}  # List of all the atomic basis vector separated by type 
    species:: Vector{String}  # Species of the atoms 
    energyFP:: Float64  # First principles energy
    modelEnergy:: Float64 # Model energy
    order::Int64 # binary, ternary, etc.
    ljvals:: Array{Float64,3}
end


struct DataSet
    crystals::Vector{Crystal}
    stdEnergy:: Float64
    meanEnergy::Float64
    offset::Float64
    nData::Int
end


#struct NS
#    L :: Float64  # Size of simulations cell (cube)
#    K :: Int64  # Number of configurations in simulation
#    W :: Int64  # Length of random walk.
#    V_ϵ :: Float64 # Convergence threshold
#    N :: Vector{Int64}  # Number of atoms per configuration
#    configs :: Vector{Crystal}  # Keeps track of all of the configurations in the simulation
##    E :: Vector{Float64}  # Energy of each configuration    
#end


# Metropolis Hasting algorithm
#struct metrop
#    nDraws:: Int
#    nBurnIn:: Int
##    μ_draws:: Array{Float64,4}
##    σ_draws:: Vector{Float64}
##    candSig_μ:: Array{Float64,3}
##    candSig_σ:: Float64
##    μ_guess:: Array{Float64,3}
##    σ_guess:: Float64
##    μ_accept:: Array{Float64,3}
##    σ_accept:: Vector{Float64}
#    proposal:: Function
#    logpost:: Function
#end
#
# General Enumeration type
struct Enum
    title:: String
    bulk:: Bool
    pLV:: Matrix{Float64}
    nD:: Int64
    dVecs:: Vector{Vector{Float64}}
    k :: Int64
    eps:: Float64

end

# Enumerated representation of a crystal
struct EnumStruct
    strN:: Int64
    hnfN::Int64
    hnf_degen:: Int64  
    lab_degen:: Int64
    tot_degen:: Int64
    sizeN:: Int64
    n:: Int64
    pgOps:: Int64
    SNF:: Diagonal{Int64,Vector{Int64}}
    HNF:: LowerTriangular{Int64, Matrix{Int64}}
    L:: Matrix{Int64}
    labeling::String
    arrows::String

end

# Molecular Dynamics Simulation
struct MD
    nParticles::Int64
    boxSize::Float64
    positions:: Array{SVector{2,Float64},1}
    velocities:: Array{SVector{2,Float64},1}
end


# The Lennard-Jones potential
struct LJ <: energyModel
    order:: Int
    params:: Array{Float64,3}
    cutoff:: Float64
end

struct LJ_metrop
    nDraws:: Int
    nBurnIn:: Int
    params_draws:: Array{Float64,4}
    σ_draws:: Vector{Float64}
    candSig_params:: Array{Float64,3}
    candSig_σ:: Float64
    params_guess:: Array{Float64,3}
    σ_guess:: Float64
    params_accept:: Array{Float64,3}
    σ_accept:: Vector{Float64}
    proposal:: Function
    logpost:: Function
    μ_Priors:: Array{Distribution,3}
    σ_Prior:: Distribution
end

# The Stillinger-Weber potential
struct SW <: energyModel
    order:: Int
    pairParams:: Array{Float64,3}
    tripParams:: Array{Float64,4}
    cutoff:: Float64
end


struct SW_metrop
    nDraws:: Int
    nBurnIn:: Int
    pairParams_draws:: Array{Float64,4}
    tripletParams_draws:: Array{Float64,4}
    σ_draws:: Vector{Float64}
    candSig_pairParams:: Array{Float64,3}
    candSig_tripletParams:: Array{Float64,3}
    candSig_σ:: Float64
    pairParams_guess:: Array{Float64,3}
    tripletParams_guess:: Array{Float64,3}
    σ_guess:: Float64
    pairParams_accept:: Array{Float64,3}
    tripletParams_accept:: Array{Float64,3}
    σ_accept:: Vector{Float64}
    proposal:: Function
    logpost:: Function
end
# Nested Sampling
struct NS
    K:: Int
    Kr:: Int
    L:: Int
    eps:: Float64
    method:: String
    configs::Vector{Crystal}
end

# Materials model
#struct model{T<:energyModel}
#    energyModel:: T
#    fitting:: metrop
#    trainingData::DataSet
#    holdoutData::DataSet
#
#end

#end