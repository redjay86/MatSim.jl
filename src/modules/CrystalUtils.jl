module CrystalUtils

using MatTypes:Crystal,DataSet,Enum
using StaticArrays
using LinearAlgebra
using StatsBase
using VASP

currentDir = abspath(@__DIR__)
libPath = relpath(joinpath(currentDir,"../libraries"))

include(abspath(libPath,"utils.jl"))
include(abspath(libPath,"structuresin.jl"))
include(abspath(libPath,"enumeration.jl"))

end # module CrystalUtils
