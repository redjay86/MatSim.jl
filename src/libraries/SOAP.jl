module SOAP


using LegendrePolynomials
using StaticArrays
using LinearAlgebra
using ase


struct descriptor
 #   crystal::Crystal.config
    n_max::Integer  # Max number of radial basis functions.
    l_max::Integer # Max angular momentum.
    r_cutoff::Float64  # Cutoff radius for the radial basis functions.
    cutoff_soft:: Float64
    cutoff_hard:: Float64
    σ::Float64  # Standard deviation for the Gaussians placed at each atom.
    W::Matrix{Float64}  # Matrix of coefficients for radial basis
#    p_vector::Vector{Float64}

end

function initialize(n_max,l_max,r_cutoff,σ)
    S = zeros(n_max,n_max)
    for α = 1:n_max,β=1:n_max
        S[α,β] = sqrt((5 + 2α) * (5 + 2β))/(5 + α + β)
        #S[α,β] = (2 * r_cutoff^(7 + α + β)) /( (5 + α + β) * (6 + α + β) * (7 + α + β))
        
    end
    W =  sqrt(LinearAlgebra.inv(S))
    return descriptor(n_max,l_max,r_cutoff,σ, W)

end

function calculate(crystal,soap)
    counter = 1
     
    p_vector = zeros(sum(0:soap.n_max) * (soap.l_max + 1)*sum(0:crystal.order))
#    for (iType, atoms)  in enumerate(crystal.positions), atom in atoms, (jType,atoms) in enumerate(crystal.positions), l = 0:soap.l_max,n=1:soap.n_max,np=1:soap.n_max
        for l = 0:soap.l_max,n=1:soap.n_max,np=n:soap.n_max
        for m  = -l:l
            c = c_nlm(crystal,crystal.positions[1],1,n,l,m,soap)  # Calculate the coefficients for the spherical harmonics.
            cp = c_nlm(crystal,crystal.positions[1],1,np,l,m,soap)  # Calculate the coefficients for the spherical harmonics.
            #println(c, cp)
            p_vector[counter] += c * conj(cp)  # Add the coefficients to the p_vector.
    
        end
        p_vector[counter] *= pi * sqrt(8/(2l + 1))
        counter += 1
    end

    return p_vector
end


function N_n(n,soap)
    return sqrt(soap.cutoff_hard/(2 * n + 5))
end

function b_n(n,r_j,soap)
    if n == -3
        return 0.0
    else
        
        if r < soap.cutoff_soft
            return -(1.0 - r/soap.cutoff_hard )^(n + 1) * soap.sigma^2/soap.cutoff_hard/N_n(n) * exp(-(r - r_j)^2 /(2 * soap.sigma^2))
                   + N_n(n - 1)/N_n(n) *(1 - r_j/soap.cutoff_hard) * b_n(n-2,r,soap)
                   + N_n(n-2)/N_n(n) * (n + 1) * soap.sigma^2/soap.cutoff_hard^2  * b_n(n-2,r,soap)
        
        elseif r < soap.cutoff_hard

        else


        end
        
        b_n(n - 1,r,soap) + b_n(n-2,r,soap)


    end



end

#Evaluate the spherical harmonics.
function Ylm_complex(l,m,θ,ϕ)
    normalizing = sqrt( (2l + 1)/(4π) * factorial(l - m)/factorial(l + m) )
    return normalizing * Plm(cos(θ),l,m) * exp(1im*m*ϕ)

end

function Ylm_real(l,m,θ,ϕ)
    if m < 0
        return sqrt(2) * imag(Ylm_complex(l,-m,θ,ϕ))
    elseif m == 0
        return Ylm_complex(l,0,θ,ϕ)
    else
        return sqrt(2) * real(Ylm_complex(l,m,θ,ϕ))
    end
end


function getLoopBounds(crystal,soap)
#    b1 = 2π * cross(crystal.lVecs[:,2],crystal.lVecs[:,3])/ dot(cross(crystal.lVecs[:,2],crystal.lVecs[:,3]),crystal.lVecs[:,1])
#    b2 = 2π * cross(crystal.lVecs[:,3],crystal.lVecs[:,1])/ dot(cross(crystal.lVecs[:,3],crystal.lVecs[:,1]),crystal.lVecs[:,2])
#    b3 = 2π * cross(crystal.lVecs[:,1],crystal.lVecs[:,2])/ dot(cross(crystal.lVecs[:,1],crystal.lVecs[:,2]),crystal.lVecs[:,3])
    b = 2π/crystal.latpar/dot(cross(crystal.lVecs[:,2],crystal.lVecs[:,3]),crystal.lVecs[:,1]) * [cross(crystal.lVecs[:,2],crystal.lVecs[:,3]) cross(crystal.lVecs[:,3],crystal.lVecs[:,1]) cross(crystal.lVecs[:,1],crystal.lVecs[:,2]) ]
    d = 2π ./ [norm(b[1]), norm(b[2]), norm(b[3])]  # Get the lengths of the reciprocal lattice vectors.
    n = convert.(Integer,ceil.(soap.r_cutoff ./ d))  # Get the number of unit cells in each direction that are needed to cover the cutoff radius.
    return n
end


function ρ(crystal,centerAtom,atomType,r,soap,separable = false)
    # Get the coordinates of the center atom in Cartesian coordinates.
    centerAtom = ase.DirectToCartesian(crystal.latpar * SMatrix(crystal.lVecs),centerAtom)
    rho = 0.0
    addVec = zeros(3)  # This is the vector that will be added to the neighbor atom to find all periodic images.
    # r is the distance from the center atom to the evaluation point.
    loopBounds = getLoopBounds(crystal,soap)  # Get the number of unit cells in each direction that are needed to cover the cutoff radius.
#    loopBounds = SVector{3,Int64}(convert.(Int64,cld.(soap.r_cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] )))    
 #   println("--------------------")
 #   println(r)
    cell = crystal.latpar * SMatrix(crystal.lVecs)
    for neighborAtom in crystal.positions[atomType] # Loop over all atoms of type "atomType" in the crystal.
        
        # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            addVec .= (i,j,k) 
            newAtom = neighborAtom .+ addVec
            newCart = ase.DirectToCartesian(cell,newAtom)
            r_j = newCart - centerAtom
            if norm(r_j) < soap.r_cutoff + 3 * soap.σ && !isapprox(norm(r_j),0,atol = 1e-3)
               # println(exp(-1 * (d^2)/(2 * soap.σ^2)))
                if separable    
                    d = norm(r_j) - norm(r)
                    d_perp = dot(r_j, norm(r_j) .* r /norm(r))
                    rho += exp(-d^2/soap.σ^2) * exp(-dot(r_j,r_j)/soap.σ^2) * exp(d_perp/soap.σ^2) 
                else
                    d = norm(r - r_j)   # This is the distance from the neighbor atom to the observation point. The gaussian should be evaluated at this distance.
                    rho += exp(-1 * (d^2)/(2 * soap.σ^2)) # 1/sqrt(2 * pi * soap.σ^2) * 
                end
            end
        end
    end 
    return rho

end

function c_nlm(crystal,centerAtom,atomType,n,l,m,soap)

    c_nlm = 0.0
    counter = 0
    dr = 0.1
    dθ = 0.1
    dϕ = 0.1
    dV = dr * dθ * dϕ  # Volume element in spherical coordinates.
    for r = 0:dr:soap.r_cutoff, θ= 0:dθ:π, ϕ = 0:dϕ:2π
#        \eps
        rCart = sphericalToCartesian(r,θ,ϕ)  # Convert spherical coordinates to Cartesian coordinates.
        Y = Ylm_real(l,m,θ,ϕ)
        #println("eval Point",rCart,r,θ,ϕ
        rhoVal = ρ(crystal,centerAtom,atomType,rCart,soap)

        c_nlm += rhoVal * 1/(exp((r - 0.8 * soap.r_cutoff)/0.2) + 1) * Y * norm(rCart)^2 * sin(θ) * gn(norm(rCart),n,soap) * dV
#        counter += 1
#        println(counter)
#        if counter > 3
#            exit("done")
#        end
    end

    return c_nlm
end

function gn(r,n,soap)
#    S = zeros(nmax,nmax)
#    for α = 1:nmax,β=1:nmax
#        S[α,β] = sqrt((5 + 2α) * (5 + 2β))/(5 + α + β)
#    end
#    W =  sqrt(LinearAlgebra.inv(S))

    g = 0.0
    for i = 1:soap.n_max
        g += soap.W[n,i] * ϕ(r,i,soap.r_cutoff)  # Calculate the radial basis function.
    end
    return g
end


function ϕ(r,n,rcut)

    N = sqrt(rcut^(2n + 5)/(2n + 5))
    return 1/N * (rcut - r)^(n + 2)
end


function sphericalToCartesian(r,θ,ϕ)
    x = r * sin(θ) * cos(ϕ)
    y = r * sin(θ) * sin(ϕ)
    z = r * cos(θ)
    return SVector(x,y,z)

end
function bispectrum(crystal,cutoff)
    loopBounds = convert.(Int64,cld.(cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] ))
    
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.positions), centerAtom in centerAtomType 
            centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
            θ = atan(centerAtomC[2]/centerAtom[1])
            ϕ = atan(centerAtom[3]/ centerAtom)
            ljvals .+= singleAtomLJ(crystal,centerAtomC,iCenter,loopBounds,params,cutoff)    # These two loops are to loop over all possible neighbor atoms
    end
end

end
