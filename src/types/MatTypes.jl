#module matTypes

#using LinearAlgebra
#using StaticArrays


abstract type energyModel end
# Crystal Structure
mutable struct Crystal
    title::String
    latpar::Float64
    lVecs:: SMatrix{3,3,Float64,9}
    nType:: Vector{Int64} #Number of each type of atom 
    aType:: Vector{Int64} # Integers representing the type for each basis atom
    nAtoms::Int64  # Total number of atoms in the unit cell
    coordSys::Vector{String} # 'D' for Direct or 'C' for cartesian
    atomicBasis:: Vector{Vector{SVector{3,Float64}}}  # List of all the atomic basis vector separated by type 
    species:: Vector{String}  # Species of the atoms 
    energyPerAtomFP:: Float64  # First principles total or peratom energy
    fitEnergy:: Float64  # Energy that should be used in the fitting process
    formationEnergyFP:: Float64  # First-principles formation energy 
    formationEnergyModel:: Float64  # model formation energy
    energyPerAtomModel:: Float64 # Model total or per atom energy
    order::Int64 # binary, ternary, etc.
    r6::Array{Float64,2}
    r12::Array{Float64,2}
end


struct DataSet
    title:: String
    crystals::Vector{Crystal}
    stdEnergy:: Float64
    meanEnergy::Float64
    offset::Float64
    nData::Int
    standardized::Bool
    fitTo::String
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

struct Pair
    atomOne :: Vector{Float64}
    atomTwo :: Vector{Float64}
    r :: Float64
    types :: Vector{Int64}
end
struct Triplet
    centerAtom:: Vector{Float64}
    atomTwo:: Vector{Float64}
    atomThree:: Vector{Float64}
    r1:: Float64
    r2:: Float64
    cosθ:: Float64
    types:: Vector{Int64}

end

# The Lennard-Jones potential
struct LJ{T}
    order:: Int64
  #  params:: Array{Float64,3}
    cutoff:: Float64

    # Replace params above.
    σ:: T#UpperTriangular{Float64, Matrix{Float64}}
    ϵ:: T#UpperTriangular{Float64, Matrix{Float64}}
    # Keep track of the standardization parameters so I can calculate the energy accurately.
    stdEnergy:: Float64  
    meanEnergy:: Float64
    offset:: Float64
    fitTo::String

end


struct LJ_metrop{D<:Distribution{Univariate,Continuous}}
    nDraws:: Int64
    nBurnIn:: Int64

    # Acceptance rates for all of the parameters
    ϵ_accept:: Array{Float64, 2}
    σ_accept:: Array{Float64, 2}
    std_accept:: Vector{Float64}
    
    # Prior distributions
    ϵ_Priors:: Array{D,2}#Vector{Distribution}
    σ_Priors:: Array{D,2}
    std_Prior:: D

    # Sigmas on proposal distributions
    ϵ_candSigs:: Array{Float64, 2}
    σ_candSigs:: Array{Float64, 2}
    std_candSig:: Float64

    # Initial guess
    std:: Float64


    # Proposal distribution
#    proposal:: F

    # Posterior
#    logpost:: G
end

# The Stillinger-Weber potential
struct SW <: energyModel
    order:: Int
    A::Array{Float64, 3}
    B::Array{Float64, 3}
    p::Array{Float64, 3}
    q::Array{Float64, 3}
    δ::Array{Float64, 3}
    a::Array{Float64, 3}
    b::Array{Float64, 3}
    λ::Array{Float64, 3}
    γ::Array{Float64,3}
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