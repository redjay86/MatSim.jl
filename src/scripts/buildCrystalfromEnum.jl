#!/usr/bin/env julia

using MatSim
structNum = try
    structNum = ARGS[1]
catch y
    error("You need to specify which structure you want to build")
end
println("Building POSCAR for structure: $structNum")

currentDir = pwd()
path = joinpath(currentDir,"struct_enum.out")
if !isfile(path)
    error("struct_enum.out not found in current directory")
end
#cd(currentDir)
enum=MatSim.read_Enum_header(path)  # Read the enum.out header
println("here")
structNum = 35
string(structNum)
enumStruct = MatSim.read_struct_from_enum(path,structNum)
print(enumStruct)
thisCrystal = MatSim.Crystal(enum,enumStruct,["Ag","Pt"])
MatSim.writePOSCAR(thisCrystal,joinpath(currentDir,"vasp." * string(structNum)))
MatSim.totalEnergy(thisCrystal,LJ_average)
display(LJ_average)