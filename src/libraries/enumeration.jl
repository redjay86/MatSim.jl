module enumeration

using LinearAlgebra:Diagonal,LowerTriangular
using Printf
#using Crystal
#using LennardJones

struct Enum
    title:: String
    bulk:: Bool
    pLV:: Matrix{Float64}
    nD:: Int64
    dVecs:: Vector{Vector{Float64}}
    k :: Int64
    eps:: Float64

end

# Enumerated representation of a crystal
struct EnumStruct
    strN:: Int64
    hnfN::Int64
    hnf_degen:: Int64  
    lab_degen:: Int64
    tot_degen:: Int64
    sizeN:: Int64
    n:: Int64
    pgOps:: Int64
    SNF:: Diagonal{Int64,Vector{Int64}}
    HNF:: LowerTriangular{Int64, Matrix{Int64}}
    L:: Matrix{Int64}
    labeling::String
    arrows::String

end




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

end


