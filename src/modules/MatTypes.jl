module MatTypes

using LinearAlgebra
using StaticArrays

currentDir = abspath(@__DIR__)
libPath = relpath(joinpath(currentDir,"../types"))

include(abspath(libPath,"MatTypes.jl"))

end # module MatTypes
