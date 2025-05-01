

cd("/Users/legoses/OneDrive - BYU-Idaho/codes/MatSim/src/libraries/")
push!(LOAD_PATH,pwd())

using PlotsMH
using LennardJones
using nSampling

using YAML
cDir = @__DIR__
cd(cDir)
# Initialize the model
resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")

# Plot the results
# Build model with averages of the draws for fit parameters.
LJ_average = PlotsMH.LJAverages(resultsPath)
LennardJones.totalEnergy(LJ_average,["Ag","Pt"])
input = YAML.load_file(joinpath(cDir,"NS.yml"))
species = ["Ag", "Pt"]
myNS = nSampling.initialize(input["params"],species,LJ_average)# Initialize the simulation...
nSampling.GMC(myNS.configs[55],500,LJ_average)
nSampling.simulate(myNS,LJ_average)
println([i.energyPerAtomModel for i in myNS.configs])
display(myNS.configs[45].coordSys)
isplay(dset.crystals[55].atomicBasis)
display(dset.crystals[55].coordSys)
display(LJ.meanEnergy)
display(model)
model = MatSim.LJ(model.order,model.cutoff, model.σ, model.ϵ,trainingSet.stdEnergy,trainingSet.meanEnergy, trainingSet.offset)


vaspDir = joinpath(cDir,"training_set")
sDir = joinpath(cDir,"structures.in")

MatSim.readVaspFolders(vaspDir,sDir,poscar = "POSCAR",energy = "total")

using YAML
using StaticArrays
using Plots
using LinearAlgebra
epsilon = 0.01
testCrystal = dset.crystals[3]
testCrystal.r6 .= 0.0
testCrystal.r12 .= 0.0
MatSim.DirectToCartesian!(testCrystal)
print(testCrystal.atomicBasis[2])
testCrystal.atomicBasis[2][1] -= SVector{3,Float64}(50 * epsilon,0,0)
loopBounds = SVector{3,Int64}(convert.(Int64,cld.(LJ_average.cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(testCrystal.latpar * testCrystal.lVecs)) )))

energies = []
forcesOne = []
forcesTwo = []
println(testCrystal.nType)
for i = 1:200
    energy = MatSim.totalEnergy(testCrystal,LJ_average)
    forceOne = MatSim.singleAtomForce(LJ_average,testCrystal,SVector(2,1),loopBounds)
    forceTwo = MatSim.gradientForce(LJ_average,testCrystal,[2,1],loopBounds)
    MatSim.DirectToCartesian!(testCrystal)
    testCrystal.atomicBasis[2][1] += SVector{3,Float64}(epsilon,0,0)
    push!(energies,energy)
    push!(forcesOne,forceOne[1])
    push!(forcesTwo,forceTwo[1])
end
length(energies)
x = -50 *epsilon:epsilon:149*epsilon
length(collect(x))
plot(x,energies)
plot(x,forcesOne)
plot(x,forcesTwo)

display(energies)


c1,c2 = MatSim.fccPures(["Ag","Pt"])
MatSim.writePOSCAR(c1,joinpath(cDir,"POSCAR.test"))
MatSim.writePOSCAR(c2,joinpath(cDir,"POSCAR2.test"))



settings = YAML.load_file(joinpath(cDir,"VASP.yml"))
kp = settings["KPOINTS"]  # KPOINTs settings
dataPath = settings["path"] # Path where the data will be found
enumSettings = settings["ENUM"] # Enum settings (indicates which structures should be built)
#enum=enumeration.read_Enum_header(enumSettings["file"])  # Read the enum.out header
potcar = settings["POTCAR"]  # POTCAR settings

# Find out how many total structures are in struct_enum.out
headerLength = length(split(readuntil(joinpath(cDir,"struct_enum.out"), "start"), "\n"))
totalStructs = countlines(joinpath(cDir,"struct_enum.out")) - headerLength

# Build the list of structures to build folders for.
if enumSettings["structs"] == "random"
    structs = sample(1:totalStructs,enumSettings["nstructs"],replace=false)
elseif enumSettings["structs"] == "sequence"
    seq = parse.(Int64,split(enumSettings["nstructs"],"-"))
    structs = seq[1]:seq[2]
end

using vaspUtils
# Build each folder.
for i in structs
    println("here: ",i)
    path = joinpath(dataPath,"str" * string(i))
    mkpath(path)
#    vaspUtils.writePOTCAR(path,settings["POTCAR"])  # Write the POTCAR file
    vaspUtils.writeINCAR(path,settings["INCAR"])   # Write the INCAR file
    #enumStruct = enumeration.read_struct_from_enum(enumSettings["file"],i)
    crystal = Crystal.fromEnum(enumSettings["file"],i,potcar["species"])
    vaspUtils.writePOSCAR(crystal,joinpath(path,"POSCAR"),)  # Write the POSCAR file.
    vaspUtils.writeKPOINTS(path,kp)
end


hull = DataSets.getConvexHull(joinpath(cDir,"gss.out"))
dPoints = [parse.(Float64,[split(x)[2],split(x)[5]]) for x in eachline(joinpath(cDir,"gss.out"))]
dPoints[:,1] 
plot([x[1] for x in dPoints],[x[2] for x in dPoints],markershape = :plus,seriestype = :scatter, label = "Convex Hull Points")
plot!(hull[:,1],hull[:,2], markershape = :circle, linecolor = :black,label = "Convex Hull")


using StaticArrays
using LinearAlgebra
SOAP.ρ(dset.crystals[152],dset.crystals[152].atomicBasis[1][3],2,SVector(0.5,0.25,0.75),20)
SOAP.c_nlm(dset.crystals[152],dset.crystals[152].atomicBasis[1][3],2,2,2,2,4)
SOAP.initialize(4,5,15)
SOAP.calculate(dset.crystals[152],SOAP.initialize(4,5,15))
display(dset.crystals[152].atomicBasis[1][1])