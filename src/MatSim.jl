module MatSim


currentDir = abspath(@__DIR__)
libPath = relpath(joinpath(currentDir,"../libraries"))
typePath = relpath(joinpath(currentDir,"../types"))

include(abspath(typePath,"MatTypes.jl"))
include(abspath(libPath,"utils.jl"))
include(abspath(libPath,"vaspUtils.jl"))
include(abspath(libPath,"metrop.jl"))
include(abspath(libPath,"enumeration.jl"))
include(abspath(libPath,"structuresin.jl"))
include(abspath(libPath,"LennardJones.jl"))


end # module MatSim
