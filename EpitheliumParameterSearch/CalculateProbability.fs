module CalculateProbability

//Abstracting the details of how everything is calculated
type technique = Simulation | Analysis

let analyticalScan allParameters = 
    let completeSet = Types.createParameterSet allParameters
    
    //Get numbers
    let probS = Array.Parallel.map (fun input -> AnalyticalCloneSizeDistribution.probabilityCloneSurvival input) completeSet
    let probN = Array.Parallel.map (fun input -> AnalyticalCloneSizeDistribution.probabilityCloneSizes input (List.init allParameters.maxN (fun i -> i+1)) allParameters.maxN) completeSet

    //Restructure results and put into a record
    AnalyticalCloneSizeDistribution.restructureParameterSet allParameters probS probN

let simulationScan allParameters =
    let completeSet = Types.createParameterSet allParameters
    let probabilityDistributions = Array.Parallel.map (fun input -> AnalyticalCloneSizeDistribution.probabilityCloneSurvival input) completeSet
    ()


let parameterSearch (input : Types.parameterSearch) approach =
    match (approach,input.deltaRange) with 
    | (Simulation,_) -> failwith "Not implemented yet"
    | (Analysis,Types.delta.Range(_)) -> failwith "Cannot generate an analytical result with non-zero delta"
    | (Analysis,Types.delta.Zero) -> analyticalScan input

