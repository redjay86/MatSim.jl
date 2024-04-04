function read_struct_from_enum(file,strN)

    f = open(file,"r")
    lines = readlines(f)
    close(f)

    title = rstrip(split(lines[1],"\n")[1])
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
    
    lineNum = strN + 15
    hnfN = parse(Int,split(lines[lineNum])[2])
    hnf_degen = parse(Int,split(lines[lineNum])[3])
    label_degen = parse(Int,split(lines[lineNum])[4])
    total_degen = parse(Int,split(lines[lineNum])[5])
    sizeN = parse(Int,split(lines[lineNum])[6])
    n = parse(Int,split(lines[lineNum])[7])
    pgOps = parse(Int,split(lines[lineNum])[8])
    SNF = Diagonal(parse.(Int,split(lines[lineNum])[9:11]))
    a = parse(Int,split(lines[lineNum])[12])
    b = parse(Int,split(lines[lineNum])[13])
    c = parse(Int,split(lines[lineNum])[14])
    d = parse(Int,split(lines[lineNum])[15])
    e = parse(Int,split(lines[lineNum])[16])
    f = parse(Int,split(lines[lineNum])[17])
    HNF = LowerTriangular([a 0 0
                           b c 0 
                           d e f])

    l = parse.(Int,split(lines[lineNum])[18:26])
    lTransform = hcat([l[i:i+2] for i=1:3:7]...)'
    labeling = split(lines[lineNum])[27]
    return Enum(title,bulk,pLattice,nD,dVecs,k,eps,strN,hnfN,hnf_degen,label_degen,total_degen,sizeN,n,pgOps,SNF,HNF,lTransform,labeling)
end