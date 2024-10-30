function read_Enum_header(file)
#    f = open(file,"r")
#    lines = readlines(f)
#    close(f)
    lines = String[]
    for line in eachline(file)
        push!(lines,string(line))
        if length(lines) >= 15
            break
        end
    end
    title = string(rstrip(split(lines[1],"\n")[1]))
    if lowercase(rstrip(split(lines[2],"\n")[1])) == "bulk"
        bulk = true
    else
        bulk = false
    end 
    bulkorsurf = rstrip(split(lines[2],"\n")[1])
    pLattice = hcat([parse.(Float64,split(split(i,"#")[1])) for i in lines[3:5]]...)  # Store vectors by columns
    nD = parse(Int,rstrip(split(lines[6])[1]))
    dVecs = [parse.(Float64,split(split(i,"#")[1])) for i in lines[7:7+nD-1]]
    k = parse(Int,rstrip(split(lines[7 + nD],"-")[1]))
    eps = parse(Float64,rstrip(split(lines[9 + nD])[1]))

    return Enum(title,bulk,pLattice,nD,dVecs,k,eps)

end


function gss(file,model,meanEnergy, stdEnergy,offset;readend = 100)
    enum=MatSim.read_Enum_header(file)
    cDir = pwd()

    io = open(joinpath(cDir,"gss.out"),"w")
    for (idx,line) in enumerate(eachline(file))
        println(idx-15)
        if idx - 15 > readend
            break
        end
        if idx < 16 
            continue
        end



        hnfN = parse(Int,split(line)[2])
        hnf_degen = parse(Int,split(line)[3])
        label_degen = parse(Int,split(line)[4])
        total_degen = parse(Int,split(line)[5])
        sizeN = parse(Int,split(line)[6])
        n = parse(Int,split(line)[7])
        pgOps = parse(Int,split(line)[8])
        SNF = Diagonal(parse.(Int,split(line)[9:11]))
        a = parse(Int,split(line)[12])
        b = parse(Int,split(line)[13])
        c = parse(Int,split(line)[14])
        d = parse(Int,split(line)[15])
        e = parse(Int,split(line)[16])
        f = parse(Int,split(line)[17])
        HNF = LowerTriangular([a 0 0
                               b c 0 
                               d e f])

        l = parse.(Int,split(line)[18:26])
        lTransform = hcat([l[i:i+2] for i=1:3:7]...)'
        labeling = split(line)[27]
        arrows = try split(line)[28] catch y repeat("0",length(labeling)) end
        eStruct =  EnumStruct(idx,hnfN,hnf_degen,label_degen,total_degen,sizeN,n,pgOps,SNF,HNF,lTransform,labeling,arrows)
        crystal = MatSim.Crystal(enum,eStruct,["Ag","Pd"],mink=false)
        energy = (MatSim.totalEnergy(crystal,model)+ offset) * stdEnergy + meanEnergy
        conc = [n/crystal.nAtoms for n in crystal.nType]
        strN = idx - 15
        printString = @sprintf "%5d  %8.4f %8.4f %8.4f\n" strN conc[1] conc[2] energy 
        write(io,printString)
    end
    close(io)
end
function read_struct_from_enum(file,strN)
    keepLine = 0
    for (idx,line) in enumerate(eachline(file))
        if idx > 15 && startswith(lstrip(line),string(strN))
            keepLine = line
            break
        end
    end
    if keepLine == 0
        error("Didn't find the specified line")
        return
    end
#    f = open(file,"r")
#    lines = readlines(f)
#    close(f)

#    title = rstrip(split(lines[1],"\n")[1])
#    if lowercase(rstrip(split(lines[2],"\n")[1])) == "bulk"
#        bulk = true
#    else
#        bulk = false
#    end 
#    bulkorsurf = rstrip(split(lines[2],"\n")[1])
#    pLattice = hcat([parse.(Float64,split(split(i,"#")[1])) for i in lines[3:5]]...)  # Store vectors by columns
#    nD = parse(Int,rstrip(split(lines[6])[1]))
#    dVecs = [parse.(Float64,split(split(i,"#")[1])) for i in lines[7:7+nD-1]]
#    k = parse(Int,rstrip(split(lines[7 + nD],"-")[1]))
#    eps = parse(Float64,rstrip(split(lines[9 + nD])[1]))
    
    lineNum = strN + 15
    hnfN = parse(Int,split(keepLine)[2])
    hnf_degen = parse(Int,split(keepLine)[3])
    label_degen = parse(Int,split(keepLine)[4])
    total_degen = parse(Int,split(keepLine)[5])
    sizeN = parse(Int,split(keepLine)[6])
    n = parse(Int,split(keepLine)[7])
    pgOps = parse(Int,split(keepLine)[8])
    SNF = Diagonal(parse.(Int,split(keepLine)[9:11]))
    a = parse(Int,split(keepLine)[12])
    b = parse(Int,split(keepLine)[13])
    c = parse(Int,split(keepLine)[14])
    d = parse(Int,split(keepLine)[15])
    e = parse(Int,split(keepLine)[16])
    f = parse(Int,split(keepLine)[17])
    HNF = LowerTriangular([a 0 0
                           b c 0 
                           d e f])

    l = parse.(Int,split(keepLine)[18:26])
    lTransform = hcat([l[i:i+2] for i=1:3:7]...)'
    labeling = split(keepLine)[27]
    arrows = try split(keepLine)[28] catch y repeat("0",length(labeling)) end
    return EnumStruct(strN,hnfN,hnf_degen,label_degen,total_degen,sizeN,n,pgOps,SNF,HNF,lTransform,labeling,arrows)
end

