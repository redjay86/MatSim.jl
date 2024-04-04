module VASP

#greet() = print("Hello World!")
using MatTypes:Crystal, DataSet
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using Printf


currentDir = abspath(@__DIR__)
libPath = relpath(joinpath(currentDir,"../libraries"))

include(abspath(libPath,"vaspUtils.jl"))

end # module VASP
