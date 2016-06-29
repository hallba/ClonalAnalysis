#load "RandomNumbers.fs"
#load "Gamma.fs"
#load "Kummer.fs"
#load "Whittaker.fs"
#load "Parallel.fs"
#load "Types.fs"
#load "IO.fs"
#load "HistoneDilution.fs"
#load "SupraBasalAnalysis.fs"
#load "SimulationCloneSizeDistribution.fs"
#load "AnalyticalCloneSizeDistribution.fs"
#load "CalculateProbability.fs"
#load "CalculateLikelihood.fs"
#load "TestData.fs"
#load "AnalyseLikelihood.fs"
#load "SupraBasalAnalysis.fs"

let sim = CalculateProbability.parameterSearch TestData.kasumip53 CalculateProbability.Simulation

let L = CalculateLikelihood.getLikelihood TestData.kasumip53Data sim

open SimulationCloneSizeDistribution

let c = {initClone with r=0.25;rho=0.5;reporting=Regular{timeLimit=4.<Types.week>;frequency=1.<Types.week>;lastReport=0.<Types.week>}}
let c2 = {c with lambda=c.lambda/3.}

let la = HistoneDilution.dilutionSinglePopulation 1000000 c |> List.map (fun i -> List.map (fun j -> log(j)) i )
let lb = HistoneDilution.dilutionTwoPopulation 660000 c 340000 c2|> List.map (fun i -> List.map (fun j -> log(j)) i )

HistoneDilution.intensityToText la "pure.dat"
HistoneDilution.intensityToText lb "mixed.dat"