module CalculateLikelihood

open MathNet.Numerics

type dataCompleteness = All | ExcludeOnes

let logFactorial i =
    //Returns log(i!) in base e
    let rec core i acc =
        if i > 0 then core (i-1) (log(float(i)) + acc) else acc/log(10.)
    core i 0.

type experimentalDataPoint = {  time: float<Types.week>
                                cloneSize: int []
                                } with 
                                member this.regularise =  ( logFactorial (Array.sum this.cloneSize) - Array.sum (Array.map (fun size -> if size = 0 then 0. else logFactorial size) this.cloneSize) )/log(10.)
                                member this.indices =   this.cloneSize
                                                        |> Array.mapi ( fun i o -> (i,o) )
                                                        |> Array.filter ( fun b -> snd(b) > 0 )
                                                        |> Array.map ( fun b -> fst(b) )
                                member this.extend n = if n > (Array.length this.cloneSize ) then 
                                                            let cloneSize' = Array.init n (fun i -> if i < (Array.length this.cloneSize ) then this.cloneSize.[i] else 0 )
                                                            {this with cloneSize=cloneSize' }
                                                       else if n < (Array.length this.cloneSize ) then failwith "Cannot curtail a set of observations"
                                                       else this
let testSystem = {  time=(11.<Types.week>/7.) ;
                    cloneSize = [| 37;13;11;6;1;4;3;1;0;1;0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|] }

//Should this be a part of the simulation code? We do not presently use it
let estimateZeroP p =
    //For zero probabilities, estimate the value based on the previous two
    //Original code does this only if the last number of clones with observations was <7
    //It then calculated the ratio of a short series of values at indices i0*.65 to i0*.75 (where i0 is the first zero)
    //Routine
    //0) Find the last index before zeros at i0
    //1) Find the series bounds i1 and i2
    //2) Find the average ratio of values in this series to the previous (wt)
    //3) Calculate future values (i.e. @ >=i) based on P(i0) * wt^(i-i0)
    let i0 = ( Array.findIndex (fun i -> i=0.) p ) - 1
    ignore (if i0 <= 7 then failwith "Hole in previous implementation" else () )
    let i1 = int(round(float(i0)*0.65))
    let i2 = int(round(float(i0)*0.75))
    let wt = Array.map2 (fun a b -> a - b)  ( Array.init (i2-i1) (fun i -> p.[i1+i]) ) ( Array.init (i2-i1) (fun i -> p.[i1+i-1]) )
             |> Array.fold (fun acc diff -> acc+diff ) 0.
             |> (fun sum -> sum/float(i2-i1))  
    Array.mapi (fun id prob -> if id < i0 then prob else p.[i0]**float((id-i0))  ) 

let logLikelihood prob obs =
    Array.map2  (fun p o -> log(p)*float(o) ) prob obs.cloneSize
    |> Array.fold (fun acc L -> L+acc ) obs.regularise

let normaliseTimePointsForSurvival cloneSizes survival =
    Array.map2 (fun nt st -> Array.map (fun nti -> nti/st) nt) cloneSizes survival

let individualLogLikelihoodContribution (pIndividual: float [] []) (search:Types.parameterSearch) (data:experimentalDataPoint list) =
    let timeMap = Map.ofArray (Array.mapi (fun i time -> (time,i)) search.timePoints)
    let correspondingP = List.map (fun dataPoint -> pIndividual.[timeMap.[dataPoint.time]] ) data //We could easily correct everything here...
    List.map2 (fun dataPoint probabilityDist -> logLikelihood probabilityDist dataPoint) data correspondingP 
    |> List.fold (fun acc p -> acc + p)  0.

let constructLogLikelihoodMatrix data P (search:Types.parameterSearch) =
    let deltaRange =    match search.deltaRange with
                        | Types.Zero -> [|0.<Types.probability>|]
                        | Types.Range(r) -> Array.ofList r
    Types.resultsMap (fun pIndividual -> individualLogLikelihoodContribution pIndividual search data ) P

let getLikelihood data (search:Types.parameterSearch) =
    let results =   match search.results with
                    | None -> failwith "Attempting to calculate a likelihood without having calculated a probability distribution"
                    | Some(res) -> res
    //Normalise count probabilities assuming survival
    let P = Types.resultsMap2 normaliseTimePointsForSurvival results.cloneSizeMatrix results.survivalMatrix
    let data = List.map (fun (obs:experimentalDataPoint) -> obs.extend search.maxN) data
    constructLogLikelihoodMatrix data P search