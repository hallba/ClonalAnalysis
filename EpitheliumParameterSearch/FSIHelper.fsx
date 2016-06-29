#load "RandomNumbers.fs"
#load "Gamma.fs"
#load "Kummer.fs"
#load "Whittaker.fs"
#load "Parallel.fs"
#load "Types.fs"
#load "IO.fs"
#load "SimulationCloneSizeDistribution.fs"
#load "AnalyticalCloneSizeDistribution.fs"
#load "CalculateProbability.fs"
#load "CalculateLikelihood.fs"
#load "TestData.fs"
#load "AnalyseLikelihood.fs"
#load "SupraBasalAnalysis.fs"

let sim = CalculateProbability.parameterSearch TestData.kasumip53 CalculateProbability.Simulation

let L = CalculateLikelihood.getLikelihood TestData.kasumip53Data sim

let c = {initClone with r=0.25;rho=0.5;reporting=Regular{timeLimit=4.<Types.week>;frequency=1.<Types.week>;lastReport=0.<Types.week>}}