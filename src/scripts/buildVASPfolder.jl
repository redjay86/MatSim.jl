#!/usr/bin/env julia

using Revise
using MatSim
using StatsBase
using YAML

cwd = pwd()
scriptdir = @__DIR__

settings = YAML.load_file(joinpath(cwd,"VASP.yml"))
kp = settings["KPOINTS"]  # KPOINTs settings
dataPath = settings["path"] # Path where the data will be found
enumSettings = settings["ENUM"] # Enum settings (indicates which structures should be built)
enum=MatSim.read_Enum_header(enumSettings["file"])  # Read the enum.out header
potcar = settings["POTCAR"]  # POTCAR settings

# Build the list of structures to build folders for.
if enumSettings["structs"] == "random"
    structs = sample(1:1633,enumSettings["nstructs"],replace=false)
end

# Build each folder.
for i in structs
    path = joinpath(dataPath,"str" * string(i))
    mkpath(path)
    MatSim.writePOTCAR(path,settings["POTCAR"])  # Write the POTCAR file
    MatSim.writeINCAR(path,settings["INCAR"])   # Write the INCAR file
    enumStruct = MatSim.read_struct_from_enum(enumSettings["file"],i)
    crystal = MatSim.Crystal(enum,enumStruct,potcar["species"])
    MatSim.writePOSCAR(crystal,joinpath(path,"POSCAR"),)  # Write the POSCAR file.
    MatSim.writeKPOINTS(path,kp)
end

pwd()
@__DIR__
cd("..")