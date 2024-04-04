using MatSim:read_struct_from_enum,Crystal,writePOSCAR

structNum = ARGS[1]
println("Building POSCAR for structure: $structNum")

currentDir = @__DIR__
cd(currentDir)
thisEnum = read_struct_from_enum("struct_enum.out",parse(Int,structNum))
thisCrystal = Crystal(thisEnum,["Ag","Pd","H"])
writePOSCAR(thisCrystal,"vasp." * structNum)

