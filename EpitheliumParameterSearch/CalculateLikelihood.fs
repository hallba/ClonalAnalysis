﻿module CalculateLikelihood

open MathNet.Numerics

type dataCompleteness = All | ExcludeOnes

let logFactorial i =
    //Returns ln(i!)
    let rec core i acc =
        if i > 0 then core (i-1) (log(float(i)) + acc) else acc
    if i > 0 then core i 0. else 0.

type experimentalDataPoint = {  time: float<Types.week>
                                cloneSize: int []
                                } with 
                                //Returns natural logarithm unless replicating matlab code- where it returns this ln(x)/(ln(10)*ln(10)) i.e. log_10(x)/ln(10)
                                member this.regularise replicateMatlab =  
                                                    let result = ( logFactorial (Array.sum this.cloneSize) - Array.sum (Array.map (fun size -> logFactorial size) this.cloneSize) )
                                                    if replicateMatlab then result / (log(10.)*log(10.)) else result
                                member this.indices =   this.cloneSize
                                                        |> Array.mapi ( fun i o -> (i,o) )
                                                        |> Array.filter ( fun b -> snd(b) > 0 )
                                                        |> Array.map ( fun b -> fst(b) )
                                member this.extend n = if n > (Array.length this.cloneSize ) then 
                                                            let cloneSize' = Array.init n (fun i -> if i < (Array.length this.cloneSize ) then this.cloneSize.[i] else 0 )
                                                            {this with cloneSize=cloneSize' }
                                                       else if n < (Array.length this.cloneSize ) then failwith "Cannot curtail a set of observations"
                                                       else this
                                member this.excludeOnes = let eOne = Array.mapi (fun id item -> if id=0 then 0 else item ) this.cloneSize
                                                          {this with cloneSize = eOne}

let testSystem = {  time=(11.<Types.week>/7.) ;
                    cloneSize = [| 37;13;11;6;1;4;3;1;0;1;0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|] }

let logLikelihood prob obs matlabReplicate =
    Array.map2  (fun p o -> if o = 0 then 0. else log(p)*float(o) ) prob obs.cloneSize
    |> Array.fold (fun acc L -> L+acc ) (obs.regularise matlabReplicate)

let normaliseTimePointsForSurvival excludeOnes (cloneSizes: float [] []) survival =
    //Exlude the probability of single cell clones
    let survival' = if excludeOnes then Array.mapi (fun i f -> f - cloneSizes.[i].[0]) survival else survival 
    Array.map2 (fun numberPatT survivalPatT -> Array.mapi (fun index specificNumberPatT -> if excludeOnes && index = 0 then 0. else specificNumberPatT/survivalPatT) numberPatT) cloneSizes survival'

let extrapolateZeroProbabilities' excludeOnes p =
    //This function is intenteded as a exact reimplementation of the matlab code
    //Where functions can return nothing (!) I've inserted failwith() terms
    //Note that erroneous zero probabilities are an issue in the original matlab implementation:
    //            % Since numerical accuracy is limited, correct zeros in P
    //            % by making use of fact that P is essentially exponential.
    let pI = Array.mapi (fun index element -> (index,element) ) p
    if (Array.filter (fun (i,e) -> ((excludeOnes && i <> 0) || not excludeOnes) && e = 0.) pI |> Array.length) = 0 then p else
        let i0 = (Array.findIndex (fun (i,e) -> i <> 0 && e=0.) pI) - 1
        if i0 <= 7 then failwith("Too few non-zero probabilities to extrapolate safely") else
            let i1 = int(round(float(i0)*0.65))
            let i2 = int(round(float(i0)*0.75))
            let a1 = Array.init (i2-i1) (fun i -> p.[i1+i-1])
            let a2 = Array.init (i2-i1) (fun i -> p.[i1+i])
            let wt = Array.map2 (fun a b -> a/b) a2 a1
                     |> Array.sum 
                     |> (fun i -> i/(float(i2-i1)))
            
            let result = Array.map (fun (i,e) -> if i <= i0 then e else p.[i0] * (pown wt (i-i0))) pI
            let resultI = Array.mapi (fun index element -> (index,element) ) result
            let zerosRemaining = (Array.filter (fun (i,e) -> ((excludeOnes && i <> 0) || not excludeOnes) && e = 0.) resultI |> Array.length)
            if zerosRemaining <> 0 then printf "%A zero elements remain\n%A\n" zerosRemaining result; failwith("Missing extrapolation")
            result

let extrapolateZeroProbabilities p =
    //pComplete is both the indices and the probabilities, stored as an array of tuples
    let pComplete = Array.mapi (fun index element -> (index,element) ) p
    if not (Array.exists (fun (index,element) -> element=0. && index <> 0) pComplete ) then p else
        //This is to replace the function from the previous implementation (see estimateZeroP)
        //This accepts an array of floats, finds the first zero and replaces all future values on 
        //the basis of an extrapolation from *at least* the last two good points

        //The value must never become 0. as a log is ultimately applied
        //To avoid this p hits a minimum it becomes the smallest possible double- 4.940656458e-324
        let i0 = ( Array.findIndex (fun (index,element) -> element=0. && index <> 0) pComplete ) - 1

        //printf "i0 = %A\n" i0
        //printf "%A %A %A %A -> i=%A\n" a b c d i0
        let lastNonZero = if i0 > 0 then p.[i0-1] else 4.940656458e-324
        //To protect against noisy simulation data the maximum ratio is 0.95, to ensure a slow decline
        //This has been arbitrarily chosen as the worst behaving datasets are the ones with an almost 
        //flat line
        let ratio =
            if i0<=0 then 
                //If we only have one or fewer values, we can't extrapolate so just set it to 0
                //printf "Everything is zero\n"
                (fun i -> 4.940656458e-324)
            else if i0 <= 7 then 
                //printf "Ratio based on %A/%A\n" p.[i0] p.[i0-1]
                let r = (p.[i0] / p.[i0-1])
                        |> fun i -> if i > 0.95 then 0.95 else i
                (fun i -> pown r i)
            else
                //try get an average ratio
                let i1 = int(round(float(i0)*0.65))
                let i2 = int(round(float(i0)*0.75))
                //printf "i1=%A i2=%A\n" i1 i2
                let r = Array.map2 (fun a b -> a/b) ( Array.init (i2-i1) (fun i -> p.[i1+i]) ) ( Array.init (i2-i1) (fun i -> p.[i1+i-1]) )
                            |> Array.fold (fun acc diff -> acc+diff ) 0.
                            |> (fun sum -> sum/float(i2-i1))  
                //printf "Series ratio %A\n" r
                let r = (p.[i0] / p.[i0-1])
                        |> fun i -> if i > 0.95 then 0.95 else i
                (fun i -> pown r i)

        Array.mapi (fun i prob -> if i < i0 then prob else 
                                                            let currentRatio= (ratio (1+i-i0))
                                                            //printf "Ratio -> %A\n" currentRatio
                                                            let a = lastNonZero*currentRatio
                                                            if a > 0. then a else 4.940656458e-324
                                                            ) p 

let individualLogLikelihoodContribution (pIndividual: float [] []) (search:Types.parameterSearch) (data:experimentalDataPoint list) =
    let timeMap = Map.ofArray (Array.mapi (fun i time -> (time,i)) search.timePoints)
    let correspondingP = List.map (fun dataPoint -> pIndividual.[timeMap.[dataPoint.time]] ) data //We could easily correct everything here...
    List.map2 (fun dataPoint probabilityDist -> logLikelihood probabilityDist dataPoint search.matlabReplicate) data correspondingP 
    //Where is the log of the multinomial   
    |> List.fold (fun acc p -> acc + p)  0.

let getLikelihood data (search:Types.parameterSearch) =
    let results =   match search.results with
                    | None -> failwith "Attempting to calculate a likelihood without having calculated a probability distribution"
                    | Some(res) -> res
    //Normalise count probabilities assuming survival & extrapolate values for "zero" probabilities
    let P = Types.resultsMap2 (normaliseTimePointsForSurvival search.excludeOnes) results.cloneSizeMatrix results.survivalMatrix 
            |> Types.resultsMap (fun pt -> Array.map (fun p -> extrapolateZeroProbabilities' search.excludeOnes p) pt)
    let data = if search.excludeOnes then List.map (fun (point:experimentalDataPoint) -> point.excludeOnes) data else data
               |> List.map (fun (obs:experimentalDataPoint) -> obs.extend search.maxN) 
    Types.resultsMap (fun pIndividual -> individualLogLikelihoodContribution pIndividual search data ) P
