#load "Gamma.fs"
#load "Kummer.fs"
#load "Whittaker.fs"
#load "Parallel.fs"
#load "Types.fs"
#load "SimulationCloneSizeDistribution.fs"
#load "AnalyticalCloneSizeDistribution.fs"
#load "CalculateProbability.fs"
#load "CalculateLikelihood.fs"
#load "TestData.fs"

let sim = CalculateProbability.parameterSearch TestData.kasumip53 CalculateProbability.Simulation

let L = CalculateLikelihood.getLikelihood TestData.kasumip53Data sim