module MD

using Random
using LinearAlgebra
using StaticArrays
using Statistics



struct MD
    nParticles::Int64
    boxSize::Float64
    positions:: Array{SVector{2,Float64},1}
    velocities:: Array{SVector{2,Float64},1}

end

function initializeParticles(nParticles::Int64,boxSize::Float64,lVecs::Array{Float64,2},random::Bool = true):: MD
    positions = zeros(SVector{2,Float64},nParticles)
    velocities = zeros(SVector{2,Float64},nParticles)
    deltaR = 0.05
    v0=1.0

    idx = 0
    for i =-boxSize-1:boxSize+1
        for j =-boxSize-1:boxSize+1
            vec = i * lVecs[1,:] + j * lVecs[2,:] + [0.6, 0.6]
            if random
                rVec = rand(Float64,2) .* 2. .- 1.
                vec +=  deltaR * rVec/norm(rVec)
            end
            if 0 <= vec[1] <= boxSize && 0 <= vec[2] <= boxSize && idx < nParticles
                idx += 1
                positions[idx] = vec

                vVec = rand(Float64,2) .* 2. .- 1.
                velocities[idx] = v0 * vVec/norm(vVec)
                println(norm(velocities[idx]), "veloc")
            end
        end
    end
    if idx < nParticles
        println("Didn't find all of the requested particles")
        println("Requested: ", nParticles)
        println("Found: ", idx)
    end
    MD(nParticles,boxSize,positions,velocities)

end

function forceOnSingleParticle(positions::Array{SVector{2,Float64},1},particle::SVector{2,Float64},boxSize::Float64):: SVector{2,Float64}
    fVec = SVector(0,0)
    #modifiedPos = SVector{2,Float64}
    #diffVec = SVector{2,Float64}
    for i=1:size(positions,1)
        diffVec = particle - positions[i]

        if abs(diffVec[1]) > boxSize/2 && abs(diffVec[2]) > boxSize/2
            modifiedPos = positions[i] + boxSize * SVector(sign(diffVec[1]) , sign(diffVec[2]))
        elseif abs(diffVec[1]) > boxSize/2
            modifiedPos = positions[i] + boxSize * SVector(sign(diffVec[1]) , 0)
        elseif abs(diffVec[2]) > boxSize/2
            modifiedPos = positions[i] + boxSize * SVector(0 , sign(diffVec[2]))
        else
            modifiedPos = copy(positions[i])
        end
        diffVec = particle - modifiedPos
        distance = norm(diffVec)
        
        if distance > 0.5
            fVec += 24 * (2/distance^13 - 1/distance^7) * diffVec/distance
        end
    end
    return fVec

end

mymod(x::SVector,y::Number) = mod.(x,y)

undowrap(x::SVector) = x .> 1

undowrapTwo(x::SVector) = x .< -1


function increaseTemp(positionsP::Array{SVector{2,Float64},1},positionsPP::Array{SVector{2,Float64},1},R::Float64)::Array{SVector{2,Float64},1}

    positionsPP = positionsP - R * (positionsP - positionsPP)

end

function simulate(myMD::MD,dt::Float64,tMax::Float64)
    R::Float64 = 1.5
    N::Int = size(myMD.positions,1)
    forces = zeros(SVector{2,Float64},N)
    T::Float64 = 0.4
    positions = copy(myMD.positions)
    positionsP = copy(myMD.positions)
    println(norm(myMD.velocities[1]), " here")
    positionsPP = positionsP - myMD.velocities * dt
    velocities = copy(myMD.velocities)
    temp = zeros(Int(tMax/dt)+1)

    idx = 1
    mask = zeros(SVector{2,Float64},N)
    for time=0:dt:tMax
        for n=1:N
            forces[n] = forceOnSingleParticle(positions,positions[n],myMD.boxSize)
        end

        @. positions = mymod(2 * positionsP - positionsPP + forces * dt^2,myMD.boxSize)

        # When a particle crosses a periodic boundary it can really
        # mess up the velocity calculation.  These next two lines
        # are to find the positions components that have changed
        # drastically and unwrap them just for the velocity calculation.

        @. mask = undowrap(positions - positionsPP) * myMD.boxSize
        @. mask += undowrapTwo(positions - positionsPP) * -myMD.boxSize
        @. velocities = (positions - positionsPP - mask)/(2 * dt)
        temp[idx] = mean(norm.(velocities).^2 ./2)

        positionsPP .= positionsP
        positionsP .= positions
        idx += 1
        if idx % 5 == 0
            plt.plot([x[1] for x in positions],[x[2] for x in positions],".")
            plt.xlim(0,4)
            plt.ylim(0,4)
            plt.draw()
            plt.pause(1e-4)
            plt.clf()
        end
        if time % 3 == 0
            println("Increasing Temp: ", R)
            @. positionsPP = positionsP - R * (positionsP - positionsPP)
        end

    end 
    
    positions,temp
end


end

#using Plots
#using Pkg
#Pkg.add("PyPlot")
#using PyPlot
#matplotlib.use("TkAgg")
##plotly()
#dt = .001
#tMax = 10.001
#lvecs =   [[1.  0.]; [0.5  sqrt(3)/2]]
#lvecs =   [[1.  0.]; [0  1]]
#myMD = initializeParticles(16,4.,lvecs)
#@time mypositions,temp = simulate(myMD,dt,tMax)
#plot(temp)
#
##plt.plot([x[1] for x in myMD.positions],[x[2] for x in myMD.positions],".")
##plt.plot(temp)
##plt.show()
#plot(map(x -> x[1],mypositions),map(x -> x[2],mypositions),seriestype = :scatter,aspect_ratio=:equal,xlim = (0,4),ylim = (0,4))
#
#using StaticArrays
#@time check = zeros(SVector{2,Float64},200)
#display(check)
#@time check[1] = SVector(1,5)
#check[2] = SVector(3,2)
#a = mean(norm.(check).^2 ./2)
#print(a)
#
#@time new = [5, 2]
#@time newtwo = SVector(-2,3)
#@time diffVec = [sign(5) , 0]
#@time new = @view SVector(sign(newtwo[1]) , 0)
#
#a = [2,4,5]
#print(a)
#
#arctan(93/245) * 180/pi