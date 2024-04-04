latpars =Dict("H"=> 3.75,"He"=> 3.57,"Li"=> 3.49,"Be"=> 2.29,"B"=> 8.73,"C"=> 3.57,"N"=> 4.039,
              "O"=> 6.83,"Ne"=> 4.43,"Na"=> 4.23,"Mg"=> 3.21,"Al"=> 4.05,"Si"=> 5.43,"P"=> 7.17,
              "S"=> 10.47,"Cl"=> 6.24,"Ar"=> 5.26,"K"=> 5.23,"Ca"=> 5.58,"Sc"=> 3.31,"Ti"=> 2.95,
              "V"=> 3.02,"Cr"=> 2.88,"Mn"=> 8.89,"Fe"=> 2.87,"Co"=> 2.51,"Ni"=> 3.52,"Cu"=> 3.61,
              "Zn"=> 2.66,"Ga"=> 4.51,"Ge"=> 5.66,"As"=> 4.13,"Se"=> 4.36,"Br"=> 6.67,"Kr"=> 5.72,
              "Rb"=> 5.59,"Sr"=> 6.08,"Y"=> 3.65,"Zr"=> 3.23,"Nb"=> 3.3,"Mo"=> 3.15,"Tc"=> 2.74,
              "Ru"=> 2.7,"Rh"=> 3.8,"Pd"=> 3.89,"Ag"=> 4.09,"Cd"=> 2.98,"In"=> 4.59,"Sn"=> 5.82,
              "Sb"=> 4.51,"Te"=> 4.45,"I"=> 7.27,"Xe"=> 6.2,"Cs"=> 6.05,"Ba"=> 5.02,"Hf"=> 3.2,
              "Ta"=> 3.31,"W"=> 3.16,"Re"=> 2.76,"Os"=> 2.64,"Ir"=> 3.84,"Pt"=> 3.92,"Au"=> 4.08,
              "Hg"=> 2.99,"Tl"=> 3.46,"Pb"=> 4.95,"Bi"=> 4.75)

element_volume =Dict("H"=>37.2958,"He"=>32.1789,"Li"=>21.2543,"Be"=>8.49323,"B"=>7.24205,"C"=>5.68741,
                 "N"=>46.6002,"O"=>22.2802,"F"=>17.0258,"Ne"=>21.7346,"Na"=>23.2596,"Mg"=>23.3928,
                 "Al"=>16.6075,"Si"=>7.8511,"P"=>9.1459,"S"=>17.1672,"Cl"=>35.2074,"Ar"=>36.3829,
                 "K"=>71.5278,"Ca"=>43.4353,"Sc"=>25.6478,"Ti"=>18.1565,"V"=>13.7718,"Cr"=>11.9439,
                 "Mn"=>19.3207,"Fe"=>11.82,"Co"=>11.1838,"Ni"=>10.9036,"Cu"=>11.7615,"Zn"=>13.311,
                 "Ga"=>18.4496,"Ge"=>19.3638,"As"=>20.4270,"Se"=>58.6173,"Br"=>33.3170,"Kr"=>46.7873,
                 "Rb"=>87.3384,"Sr"=>56.1889,"Y"=>33.0792,"Zr"=>23.8327,"Nb"=>17.9685,"Mo"=>15.6279,
                 "Tc"=>14.5458,"Ru"=>13.9206,"Rh"=>13.718,"Pd"=>14.716,"Ag"=>17.1045,"Cd"=>18.7161,
                 "In"=>26.6861,"Sn"=>29.3238,"Sb"=>27.1733,"Te"=>62.3227,"I"=>24.3807,"Xe"=>59.582,
                 "Cs"=>110.723,"Ba"=>63.253,"Hf"=>23.1748,"Ta"=>18.1323,"W"=>15.7772,"Re"=>14.8694,
                 "Os"=>14.5485,"Ir"=>14.1558,"Pt"=>15.0591,"Au"=>16.9793,"Hg"=>27.6914,"Tl"=>29.2949,
                 "Pd"=>30.3218,"Bi"=>31.2849,"U"=>13.6389, "O"=>13.6389)#"U"=>50.84, "O"=>15.86}


Crystal(folder::String,file::String, species::Vector{String}; overwriteLatPar = false,energyFP = 0, energyPred = 0) = fromPOSCAR(folder,file,species,overwriteLatPar = overwriteLatPar, energyFP = 0, energyPred = 0)
Crystal(list::Vector{String}, species::Vector{String}; overwriteLatPar = false, energyFP = 0, energyPred = 0) = fromPOSCAR(list, species,overwriteLatPar = overwriteLatPar,energyFP = energyFP, energyPred = energyPred)
Crystal(enum::Enum,species::Vector{String}) = fromEnum(enum,species,energyFP = 0,energyPred = 0,ljvals = zeros(3,2))

function fromEnum(enum::Enum,species:: Vector{String};energyFP = 0, energyPred = 0,ljvals = zeros(3,2))
    sLV = enum.pLV * enum.HNF  # Get super lattice vectors!
    a,b,c,d,e,f = enum.HNF[1,1],enum.HNF[2,1],enum.HNF[2,2],enum.HNF[3,1],enum.HNF[3,2],enum.HNF[3,3]
    # Get atomic basis..
    nAtoms = length(enum.labeling)
    aBas = [zeros(3) for i =1:nAtoms]
    gIndex = zeros(Int64,nAtoms)
    iC = 0
    for iD = 1:enum.nD
        for z1 = 0:a-1, z2= Int((b * z1) / a):Int(c + (b * z1) / a) - 1, z3 = Int(z1 * (d - (e * b) / c) / a + (e * z2) / c):Int(f + z1 * (d - (e * b) / c) / a + (e * z2) / c) - 1
            iC += 1
            temp = enum.pLV * [z1,z2,z3]
            temp2 = temp + enum.dVecs[iD]

            aBas[iC] .= temp2

            gReal = enum.L * [z1,z2,z3]

            g = Int.(enum.L * [z1,z2,z3])
            if !(g ≈ gReal)
                error("Mapping didn't work")
            end

            # Bring g-vector back into the first cell
            g = g .% diag(enum.SNF)
            # Which character in the labeling string corresponds to this atom..
            gIndex[iC] = Int((iD - 1) * enum.SNF[1,1] * enum.SNF[2,2] * enum.SNF[3,3]
                         + g[1] * enum.SNF[2,2] * enum.SNF[3,3]
                         + g[2] * enum.SNF[3,3]
                         + g[3]) + 1
            
        end
    end
    aType = [parse(Int,enum.labeling[i]) for i in gIndex]  # Get the atomic labels in the right order (i.e. using gIndex)
    nType = [count(==(i),aType) for i in Set(aType)]

    if !all(isapprox.(sort(gIndex), gIndex,atol = 1e-7))
        println(aType)
        println(sort(aType))
        println(gIndex)
        #println("Found one")
        error("This isn't an error. I was just curious if the gIndex Vector is ever in an odd order.")
    end
    idx = sortperm(aType)  # Sort the type list, saving the indices
    aBas = getPartitions(aBas[idx],nType) #Sort the atomic basis to match and then partition them by atom type.
    aType = sort(aType)

    cellVolume = abs(cross(sLV[:,1],sLV[:,2])' * sLV[:,3])
    latpar = vegardsVolume(species,nType,cellVolume)
    final = Crystal(enum.title * " str #: " * string(enum.strN),latpar,minkowski_reduce(sLV,1e-7),nType,aType,nAtoms,["C"],aBas,["Unk" for i in nType],0,0,enum.k,zeros(Float64,3,2))
    CartesianToDirect!(final)
    return final
end

function getPartitions(atoms::Vector{Vector{Float64}},nAtoms::Vector{Int64})
    if length(atoms) != sum(nAtoms)
        println("Number of atoms doesn't match")
        return
    end
    parts = prepend!(cumsum(nAtoms),0)
    return [atoms[parts[i]+1:parts[i+1]] for i=1:length(parts)-1]
end

function buildRandom(lPar:: Float64, lVecs::Matrix{Float64},nAtoms::Vector{Int64},cutoff::Float64)

    totalAtoms = sum(nAtoms)
    order = length(nAtoms)
    aType = hcat([ [n for i=1:nAtoms[n]]' for n=1:length(nAtoms)]...)'

    atoms = [zeros(Float64,3) for i=1:totalAtoms]  

    nTotal = 0
    counter = 0
    while nTotal < totalAtoms
        counter += 1
        if counter > 5000
            println("Having trouble putting that many atoms into this simulation cell.")
            println("So far I've only place $nTotal atoms." )
            return
        end
        newAtomCart = lPar * lVecs * rand(3)

        newOK = true
        for i=1:nTotal
            if norm(newAtomCart - atoms[i]) < cutoff
                newOK = false
            end
            if !newOK
                break
            end
        end
        if newOK
            nTotal += 1
            atoms[nTotal] .= newAtomCart
        end
    end
    atomicBasis = getPartitions(atoms,nAtoms)
    return Crystal("Random Locations",lPar,lVecs,nAtoms,aType,totalAtoms,["C"], atomicBasis,0.0,0.0,order,zeros(3,2))
end


function fromPOSCAR(folder::String,file::String,species::Vector{String};overwriteLatPar = false, energyFP = 0, energyPred = 0)

    cd(folder)
    
    file = open(file, "r")
    pos = readlines(file)

    title = pos[1]
    lVecs = SMatrix{3,3}(reduce(hcat,[parse.(Float64, y) for y in [split(x) for x in pos[3:5]]])) # Read in the lattice vectors.
    if !isletter(lstrip(pos[6])[1])
        counter = 6
    else
        counter = 7
    end
    nBasis = parse.(Int,split(pos[counter]))
    coordSys = [pos[counter + 1]]
    latpar = parse(Float64,pos[2])
    aType = hcat([ [n for i=1:nBasis[n]]' for n=1:length(nBasis)]...)'
    println(pos[8:7 + sum(nBasis)])
    allBasis = [SVector{3,Float64}(parse.(Float64,split(x)[1:3])) for x in pos[(counter + 2):(counter + 1 +sum(nBasis))]] # Read all of the basis vectors
    allTypes = try
        [split(x)[end] for x in pos[(counter + 2):end]]
    catch e
        println("No atom types listed in the POSCAR")
        ["?" for x in pos[(counter + 2):end]]
    end
    atomicBasis = getPartitions(allBasis,nBasis)
    order = length(nBasis)
    nAtoms = sum(nBasis)


    # I'm not sure why I'm putting this code here. We'll always want to use
    # the lattice parameter found in the structures.in file because that
    # was the geometry used in the VASP calculation. 
    cellVolume = abs(cross(lVecs[:,1],lVecs[:,2])' * lVecs[:,3])
    if overwriteLatPar
        println("Overwriting lat par in file with Vegard's law")
        latpar = vegardsVolume(species,nBasis,cellVolume)
    else
        println("Keeping lattice parameter in file")
    end
    
    

    return Crystal(title, latpar,lVecs,nBasis,aType,nAtoms,coordSys,atomicBasis,species,energyFP,energyPred,order,zeros(3,2))  # Create new crystal object.

end




function fromPOSCAR(lines::Vector{String},species::Vector{String};overwriteLatPar = false, energyFP = NaN,energyPred = NaN)
    title = lines[1]
    lVecs = SMatrix{3,3}(reduce(hcat,[parse.(Float64, y) for y in [split(x) for x in lines[3:5]]])) # Read in the lattice vectors.
    if !isletter(lstrip(lines[6])[1])
        counter = 6
    else
        counter = 7
    end
    nBasis = parse.(Int,split(lines[counter]))
    coordSys = [lines[counter + 1]]
    latpar = parse(Float64,lines[2])
    aType = hcat([ [n for i=1:nBasis[n]]' for n=1:length(nBasis)]...)'
    allBasis = [parse.(Float64,split(x)[1:3]) for x in lines[(counter + 2):(counter + 1)+sum(nBasis)]] # Read all of the basis vectors
    allTypes = try
        [split(x)[end] for x in lines[(counter + 2):end]]
    catch e
        println("No atom types listed in the POSCAR")
        ["?" for x in pos[(counter + 2):end]]
    end
    atomicBasis = getPartitions(allBasis,nBasis)
    order = length(nBasis)
    nAtoms = sum(nBasis)

    # I'm not sure why I'm putting this code here. We'll always want to use
    # the lattice parameter found in the structures.in file because that
    # was the geometry used in the VASP calculation. 
    cellVolume = abs(cross(lVecs[:,1],lVecs[:,2])' * lVecs[:,3])
    if overwriteLatPar
        println("Overwriting lat par in file with Vegard's law")
        latpar = vegardsVolume(species,nBasis,cellVolume)
    else
        println("Keeping lattice parameter in file")
    end
    return Crystal(title, latpar,lVecs,nBasis,aType,nAtoms,coordSys,atomicBasis,species,energyFP,energyPred,order,zeros(3,2))  # Create new crystal object.

end



function DirectToCartesian!(crystal::Crystal)
    if crystal.coordSys[1] == "D"
        println("Converting to cartesian")
        crystal.atomicBasis .= [[crystal.lVecs * i for i in j] for j in crystal.atomicBasis ]
        crystal.coordSys[1] = "C"
    else
        println("Already in Cartesian coordinates")
    end
end

function CartesianToDirect!(crystal::Crystal)
    if crystal.coordSys[1] == "C"
        println("Converting to direct")
        crystal.atomicBasis .= [[ round.( ( inv(crystal.lVecs) * i) .% 1,sigdigits = 8) for i in j]  for j in crystal.atomicBasis ]

        crystal.coordSys[1] = "D"
    else
        println("Already in Cartesian coordinates")
    end
end

function DirectToCartesian(lVecs::Matrix{Float64},atom:: Vector{Float64})
    return lVecs * atom
end

function CartesianToDirect(lVecs::Matrix{Float64},atom:: Vector{Float64})
    return inv(lVecs) * atom  .% 1
end

function isEquivalentAtom(atomOne,atomTwo)  # Atoms are assumed to be in direct coordinates.

    #=If the atoms are separated by a lattice vector we can be sure that they are equivalent atoms (translationally).
       To check this, subtract the two vectors and see if the resulting vector is only composed of 
       integers. 
    =#
    return all(round.(Int64,atomOne - atomTwo) - (atomOne - atomTwo) == [0,0,0]) 


end




function Gaussian_Reduce(U,V,eps)  
    """This routine takes two vectors (in three-space) and reduces them to
    form a shortest set (Minkowski reduced). The idea is to subtract
    multiples of U from V so that the new V is as close to the origin
    as any lattice point along the line that passes through U in the
    direction of V. The process is repeated until the new vector isn't
    shorter than the other. It's pretty obvious if you do an example
    by hand. Also see 3.1 of Lecture notes in computer science, ISSN
    0302-974, ANTS - VI: algorithmic number theory, 2004, vol. 3076,
    pp. 338-357 ISBN 3-540-22156-5. Fixes made Apr 2012 GLWH (not sure
    if they made a practical difference though)
    :arg U: a vector
    :arg V: a vector
    :arg eps: finite precision tolerance
    """

#    from numpy.linalg import norm
#    from numpy import dot

    it = 0
    if norm(U) > (norm(V) - eps)
        # Make sure that the {U,V} are listed in ascending order; ||U||<||V||
        temp = copy(U)
        U = copy(V)
        V = copy(temp)  # Keep V as the longest vector
    end
    done = false
    it = 1
    while  !done
        if it > 10  # pragma: no cover
            error("gaussian_reduce_two_vectors failed to converge in 10 iterations")
        end
        R = V - Int(round((U'*V)/(U'*U) + 1e-10)) * U

        V = copy(U)  # Swap U and V (so U remains the shortest)
        U = copy(R)
        if norm(U) >= (norm(V) - eps)
            done = true
        end
        it += 1
    end
    # Make sure that the {U,V} are listed in ascending order on exit; ||U||<||V||
    temp = copy(U)
    U = copy(V)
    V = copy(temp)
    return U, V
end


function reduce_C_in_ABC(ABC,eps)
    #oldABC = deepcopy(ABC)
    A,B,C = ABC[:,1],ABC[:,2], ABC[:,3] 
    A,B = Gaussian_Reduce(A,B,eps)
    cpdAB = cross(A,B)/norm(cross(A,B))  # Find the unit vector the points perp to A-B plane

    T = C - cpdAB * (C' * cpdAB)  # Find projection of C onto A-B plane
    if !isapprox(T' * cross(A,B),0,atol = eps)  # Check to see if we did it right. Dot product should be zero now.
        error("Projection of C into A,B plane failed.")
    end

#    ABC = hcat(A, B, C)  # lVs by columns, not rows.
    LC = Int.(floor.(inv(ABC) * T .+ eps)) # Find what combinations of A,B,C will produce the projection.
    corners = [0 0 0
               1 0 0
               0 1 0
               1 1 0]  # The projection is closest to one of these corners.
    distances = [norm(T - ABC * (corners[i,:] + LC) ) for i =1:4]  # Get the distances to each of the corners, so I know which combination of ABC to subtract off.
    val,idx = findmin(distances)

    C = C - (ABC * (corners[idx,:] + LC))  # This is guaranteed to produce another lattice vector because you are subtracting multiples of other lattice vectors.

    newABC = hcat(A,B,C)

    if !all(isapprox.(inv(newABC) * ABC - Int.(round.(inv(newABC) * ABC)),0,atol = eps))
        display(inv(newABC)* ABC)
        display(Int.(round.(inv(newABC) * ABC)))
        println("This matrix should have been integers")
        error("Lattice not preserved in reduce_C_in_ABC")
    end
    return newABC

end


function minkowski_reduce(ABC,eps)
    #ABC = hcat(A,B,C)
    outABC = deepcopy(ABC)
    limit = 10
    for it=1:limit
        norms = [norm(x) for x in eachcol(outABC)]
        outABC = outABC[:,sortperm(norms)]
        outABC = reduce_C_in_ABC(outABC,1e-7)
        if norm(outABC[:,3]) >= norm(outABC[:,2]) - eps
            break
        end
    end

    if !minkowski_check(outABC,eps)
        error("Minkowski conditions not met")
    end
    if det(outABC') < 0
        save = outABC[:,2]
        outABC[:,2] .= outABC[:,3]
        outABC[:,3] .= save
    end
    return outABC
end


function minkowski_check(ABC,eps)

    b1 = ABC[:,1]
    b2 = ABC[:,2]
    b3 = ABC[:,3]

    minkowski_check = true

    if norm(b1) > norm(b2) + eps
        minkowski_check = false
        println("Minkowski_condition 1 failed: b1 > b2")
    end

    if norm(b2) > norm(b3) + eps
        minkowski_check = false
        println("Minkowski_condition 2 failed: b2 > b3")
    end

    if norm(b2) > norm(b1 + b2) + eps
        minkowski_check = false
        println("Minkowski_condition 3 failed: b2 > b1 + b2")
    end

    if norm(b2) > norm(b1-b2) + eps
        minkowski_check = false
        println("Minkowski_condition 4 failed: b2 > b1-b2")
    end

    if norm(b3) > norm(b1+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 5 failed: b3 > b1+b3")
    end

    if norm(b3) > norm(b3-b1) + eps
        minkowski_check = false
        println("Minkowski_condition 6 failed: b3 > b3-b1")
    end

    if norm(b3) > norm(b2+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 7 failed: b3 > b2+b3")
    end

    if norm(b3) > norm(b3-b2) + eps
        minkowski_check = false
        println("Minkowski_condition 8 failed: b3 > b3-b2")
    end

    if norm(b3) > norm(b1+b2+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 9 failed: b3 > b1+b2+b3")
    end

    if norm(b3) > norm(b1-b2+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 10 failed: b3 > b1-b2+b3")
    end

    if norm(b3) > norm(b1+b2-b3) + eps
        minkowski_check = false
        println("Minkowski_condition 11 failed: b3 > b1+b2-b3")
    end

    if norm(b3) > norm(b1-b2-b3) + eps
        minkowski_check = false
        println("Minkowski_condition 12 failed: b3 > b1-b2-b3")
    end

    return minkowski_check
end






function vegardsVolume(elements,atom_counts,volume)
    if nothing in [findfirst(x-> x == lowercase(i),lowercase.(keys(element_volume))) for i in elements]
        return 1.0
    end
    
    nAtoms = sum(atom_counts)
    nTypes = length(atom_counts)
    concentrations = [x/nAtoms for x in atom_counts]
    if length(elements) != nTypes
        error("Not enough elements specified")
    end

    #  The next two lines are Wiley's way.
    # Find lat pars if crystal was pure.
    purelatPars = [( element_volume[elements[x]]/(volume/nAtoms) )^(1. /3.) for x =1:nTypes]
    # Concentration-weighted average over lattice parameters
    wiley = sum([purelatPars[x] * concentrations[x] for x =1:nTypes])


    # The next two lines are my way.
    # concentration-weighted average of volumes
    avgVolume = sum([element_volume[elements[x]] *concentrations[x] for x =1:nTypes] )
    # Calculate lattice parameter
    mine = ( avgVolume/(volume/nAtoms) )^(1/3.)
#    print("My approach: {}.  Wiley's approach: {}".format(mine,wiley))
    return wiley
end
