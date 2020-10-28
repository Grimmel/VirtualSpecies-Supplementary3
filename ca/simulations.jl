include("ca.jl")
using .CellularAutomata
using DelimitedFiles
using Profile
function runSimulations()
    # Parameters
    landscapes = ("LLM","LMS","LLS","LMM","LSS")
    landscapeReplicates = 10
    simulationReplicates = 9 # 9 is the max, 10 total
    speciesList = Dict("ca_species1"=>[0.6,2],
                        "ca_species2"=>[0.6,4],
                        "ca_species3"=>[0.8,2],
                        "ca_species4"=>[0.8,4])
    for (species,params) in speciesList
        for landscape in landscapes
            for lsrep in 1:10
                for srep in 0:9
                    suit = readdlm("D:/PHDExperimentOutputs/SimLandscapes/suitability/"*landscape*string(lsrep)*"_suitability.asc",skipstart=6)
                    suit = suit ./100
                    ls_dimension = size(suit)
                    # Dispersal Parameters
                    dispersalProbability = params[1]
                    dispersalDistance = params[2]
                    n_iterations = 200
                    pa = ones(400,400)
                    pa_cart_index = CartesianIndices(pa)
                    simulationModel = OccurenceCellularAutomata(pa,pa_cart_index,suit,dispersalProbability,0,0,0,dispersalDistance)
                    for i in 1:n_iterations
                        colonise(simulationModel)
                        extinction(simulationModel)
                    end
                    open("D:/PHDExperimentOutputs/MainSims/"*species*"/Output_maps/abundance/abundance_s1_"*landscape*string(lsrep)*"_r"*string(srep)*".tif","w") do io
                        write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
                        writedlm(io,simulationModel.pa)
                    end
                end
            end
        end
    end
end
