module CalculateLikelihood

open MathNet.Numerics

type dataCompleteness = All | ExcludeOnes

let logFactorial i =
    //Returns log(i!) in base e
    let rec core i acc =
        if i > 0 then core (i-1) (log(float(i)) + acc) else acc/log(10.)
    core i 0.

type experimentalDataPoint = {  time: float<CloneSizeDistribution.week>
                                cloneSize: int []
                                } with 
                                member this.regularise =  ( logFactorial (Array.sum this.cloneSize) - Array.sum (Array.map (fun size -> if size = 0 then 0. else logFactorial size) this.cloneSize) )/log(10.)
                                member this.indices =   this.cloneSize
                                                        |> Array.mapi ( fun i o -> (i,o) )
                                                        |> Array.filter ( fun b -> snd(b) > 0 )
                                                        |> Array.map ( fun b -> fst(b) )

let testSystem = {  time=(11.<CloneSizeDistribution.week>/7.) ;
                    cloneSize = [| 37;13;11;6;1;4;3;1;0;1;0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|] }



let getLikelihood data (probabilities:CloneSizeDistribution.parameterSearch )=
    let C = Array.map (fun (dataPoint:experimentalDataPoint) -> dataPoint.regularise) data
    let ind = Array.map (fun (dataPoint:experimentalDataPoint) -> dataPoint.indices) data

    let P = Array.map2 (fun (surv:complex) size -> Array.map (fun oneSize -> oneSize/surv.r) size) probabilities.oneDimSurvMatrix probabilities.oneDimSizeMatrix

    0.