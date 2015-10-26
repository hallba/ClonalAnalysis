module CalculateProbability

//Abstracting the details of how everything is calculated
type technique = Simulation | Analysis

let analyticalScan allParameters = 
    let completeSet = Types.createParameterSet allParameters
    
    //Get numbers
    let probS = Array.Parallel.map (fun input -> (AnalyticalCloneSizeDistribution.probabilityCloneSurvival input).Real) completeSet
    let probN = Array.Parallel.map (fun input -> AnalyticalCloneSizeDistribution.probabilityCloneSizes input (List.init allParameters.maxN (fun i -> i+1)) allParameters.maxN) completeSet

    //Restructure results and put into a record
    AnalyticalCloneSizeDistribution.restructureParameterSet allParameters probS probN

let simulationScan allParameters number =
    let completeSet = Types.createParameterSet allParameters
    let probabilityDistributions = Array.Parallel.map (fun input -> input 
                                                                    |> SimulationCloneSizeDistribution.parameterSetToClone (SimulationCloneSizeDistribution.Specified(List.ofArray allParameters.timePoints))
                                                                    |> SimulationCloneSizeDistribution.cloneProbability number
                                                                    ) completeSet
    let oneDLength = (Array.length probabilityDistributions) * (Array.length allParameters.timePoints)
    let convert i = ( (i % (Array.length allParameters.timePoints)), (i / (Array.length probabilityDistributions)) )
    let probabilityDistributions1D = Array.init oneDLength (fun i -> let (t,r) = convert i
                                                                     probabilityDistributions.[r].[t].basalFraction )
    //Todo probS is 1 - probabilityDistributions.[x].[0]
    let probS = Array.Parallel.map (fun (i: float []) -> 1. - i.[0]) probabilityDistributions1D
    //Todo probN is probabilityDistributions.[x][1..]
    let probN = Array.Parallel.map (fun (i: float []) -> Array.init ((Array.length i)-1) (fun j -> i.[j]) ) probabilityDistributions1D
    AnalyticalCloneSizeDistribution.restructureParameterSet allParameters probS probN


let parameterSearch (input : Types.parameterSearch) approach =
    match (approach,input.deltaRange) with 
    | (Simulation,_) -> failwith "Not implemented yet"
    | (Analysis,Types.delta.Range(_)) -> failwith "Cannot generate an analytical result with non-zero delta"
    | (Analysis,Types.delta.Zero) -> analyticalScan input

