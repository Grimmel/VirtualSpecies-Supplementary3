module CellularAutomata
using Distributions
using DelimitedFiles
using Random
using StatsBase

export colonise,coloniseNeighbourWeight,extinction,extinctionNeighbour,OccurenceCellularAutomata


struct OccurenceCellularAutomata
    pa::Matrix{Float64}
    pa_cart_index::Matrix{CartesianIndex{2}}
    suitability::Matrix{Float64}
    dispersalProbability::Float64
    meanDispersal::Float64
end

function newPos(meanDistance::Float64,cartIdx::CartesianIndex)
    distance = rand(Exponential(meanDistance),1)
    # + 1 ensures dispersal outside of the initial cell
    distance = distance[1] + 0.75
    angle = 360.0*rand()
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + cartIdx[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + cartIdx[1]
    coord = Tuple(y::Int64,x::Int64)
    return(coord)
end
function selectProportion(pa::Array{Float64},caIndex::Array{CartesianIndex},dispersalProbability::Float64)
    nPresences = Int(sum(pa))
    total = Array{CartesianIndex}(undef, nPresences)
    counter = 1
    for idx in caIndex
        if pa[idx] === 1.0
            total[counter] = idx
            counter+=1
        end

    end
    numberDispersing = sample(total,rand(Binomial(nPresences,dispersalProbability)))
end
function colonise(ca::OccurenceCellularAutomata)
    shape = size(ca.pa)
    dCells = selectProportion(ca.pa,ca.pa_cart_index,ca.dispersalProbability)
    for i in dCells
        newXY = newPos(ca.meanDispersal,i)
        if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
            ca.pa[newXY[1],newXY[2]] = 1.0
        end
    end
end
function extinction(ca::OccurenceCellularAutomata)
    nPresences = Int(sum(ca.pa))
    total = Array{CartesianIndex}(undef, nPresences)
    for idx in ca.pa_cart_index
        if ca.pa[idx] === 1.0
            survived = rand(Bernoulli(ca.suitability[idx]),1)
            if survived[1] === false
                ca.pa[idx] = 0.0
            end
        end
    end
end
end
