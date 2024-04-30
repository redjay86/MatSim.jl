#module VASP





#function readStructuresIn(folder::String,file::String,species::Vector{String};overwriteLatPar = false)
#    cd(folder)
#    file = open(file,"r")
#    pos = readlines(file)
#
#    data = Vector{Crystal}()
#    for (idx,line) in enumerate(pos)
#
#        if occursin("#--",line)
#            nAtoms = sum([parse(Int64,x) for x in split(pos[idx + 6])])
#            startpoint = idx + 1
#            theend = idx + 7 + nAtoms
#            thisCrystal = Crystal(pos[startpoint:theend],species,overwriteLatPar = overwriteLatPar, energyFP = parse(Float64,pos[theend + 2 ]))
#            if !isnan(thisCrystal.energyFP)
#                push!(data,thisCrystal)
#            end
#        end
#    end
#    return DataSet(data)
#end

function readVaspFolders(folder::String,file::String;poscar = "CONTCAR",outcar = "OUTCAR")

    for obj in readdir(folder,join = true)

        if isdir(obj)
            println(obj)
            poscarPresent = isfile(joinpath(obj,poscar))
            outcarPresent = isfile(joinpath(obj,outcar)) 
            if poscarPresent && outcarPresent
                crystal = Crystal(obj,poscar,energyFP = getEnergy(joinpath(obj,outcar)))
                writePOSCAR(crystal,joinpath(folder,file),"a")
                open(joinpath(folder,file), "a") do f
                    write(f,"Energy:","\n")
                    writedlm(f,crystal.energyFP)
                    write(f,"#--------------------","\n")
                end
                
                
            else
                @printf("Didn't find one of the necessary files. POSCAR: %s OUTCAR: %s Skipping this dir: %s\n",poscarPresent ? "true" : "false",outcarPresent ? "true" : "false", obj)
                #print(obj)
            end
        end
    end



end


function writePOSCAR(crystal::Crystal,fName::String,openCode::String = "w")
    open(fName,openCode) do f
        write(f,crystal.title,'\n')
        writedlm(f,crystal.latpar',' ')
        writedlm(f,crystal.lVecs',' ')
        writedlm(f,crystal.nType',' ')
        write(f,crystal.coordSys[1],'\n')
        counter = 1
        for basis in crystal.atomicBasis
            for atom in basis
                writedlm(f,atom',' ')
               # seek(f,position(f) - 1)
            #    write(f,' ' * crystal.atomTypes[counter],'\n')
                counter +=1
            end
        end
    
    end
    

end

function getEnergy(filePath::String)

    open(filePath) do f
        lines = readlines(f)
        energyLine = findall([occursin("free  energy",x) for x in lines])
        return parse(Float64,split(lines[energyLine[1]])[end-1])
    end

end

function writePOTCAR(path::String,potcar::Dict)

#    settings = YAML.load_file("VASP.yml")
    potcarsroot = potcar["path"]
    potcarspecies = reverse(sort(potcar["species"]))
    n = length(potcarspecies)
    dirs = [potcarsroot * "/" * x for x in potcarspecies]
        
    catCommand = `cat $dirs`
    outpath = joinpath(path,"POTCAR")
    run(pipeline(catCommand, stdout = outpath))

end


function writeINCAR(writepath::String,incar::Dict)

    open(joinpath(writepath,"INCAR"), "w") do f
        for (i,d) in incar
            write(f,i * "=" * string(d) * "\n")
        end
    end
end

function writeKPOINTS(writepath::String,kpoints::Dict)

    # generate the input file.

    open(joinpath(writepath,"KPGEN"),"w") do f
        for (key,val) in kpoints["settings"][1]
            write(f,key * "=" * string(val)*"\n")
        end
        write(f,"")
    end
#    potcarsroot = potcar["path"]
#    potcarspecies = reverse(sort(potcar["species"]))
#    n = length(potcarspecies)
#    dirs = [potcarsroot * "/" * x for x in potcarspecies]
#        
    currentDir = pwd()
    cd(writepath)
    kpCommand = `kpoints.x`
    run(kpCommand)
    #wait()
    cd(currentDir)
#    outpath = joinpath(path,"POTCAR")
#    run(pipeline(catCommand, stdout = outpath))

end

#end

#myElements = ["Ni", "Pd"]
#
#atoms = [5,8]
#volume = 105
#vegardsVolume(myElements,atoms,volume)