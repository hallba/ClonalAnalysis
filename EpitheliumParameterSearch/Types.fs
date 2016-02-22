module Types

open MathNet.Numerics
open System

[<Measure>]
type cell

[<Measure>]
type week

[<Measure>]
type probability

type delta = Zero | Range of float<probability> list

type scanResults       = {  cloneSizeMatrix:  float [] [] [] [] [] []  // delta x lambda x rho x r x time x n
                            survivalMatrix: float [] [] [] [] []  //delta x lambda x rho x r x time
                            oneDimSizeMatrix: float [] [] 
                            oneDimSurvMatrix: float [] 
                            }

[<Serializable>]
type parameterSearch = {    rhoRange:       float array
                            rRange:         float<probability> array
                            lambdaRange:    float<cell/week> array
                            timePoints:     float<week> array
                            maxN:           int
                            deltaRange:     delta
                            matlabReplicate:bool
                            supraBasalFit:  bool
                            results:        scanResults option
                            excludeOnes:    bool
                            }

type parameterSet = {   time: float<week>;
                        rho: float; 
                        r: float<probability>;
                        delta: float<probability>
                        lambda: float<cell/week>
                        maxN:   int
                        } with 
                        member this.migration = this.rho/(1.-this.rho) * this.lambda
                        member this.correctPathologicalPoint = { this with rho = this.rho + 0.0001 ; r = this.r+0.00001<probability> }
                        override this.ToString() = sprintf "Rho: %A r:%A Time:%A Lambda:%A" this.rho this.r this.time this.lambda

let testSystem = {time=1.<week>; rho=0.85; r=0.15<probability>; delta=0.<probability>; lambda=2.<cell/week>; maxN=10}

let createParameterSet (input : parameterSearch) = 
    match input.deltaRange with
    | Zero ->
                [ for lambda in input.lambdaRange do for rho in input.rhoRange do for r in input.rRange do for t in input.timePoints do yield {testSystem with lambda=lambda;time=t;r=r;rho=rho;maxN=input.maxN} ]
                |> Array.ofList
    | Range(d) ->
                [ for delta in d do for lambda in input.lambdaRange do for rho in input.rhoRange do for r in input.rRange do for t in input.timePoints do yield {testSystem with delta=delta;lambda=lambda;time=t;r=r;rho=rho;maxN=input.maxN} ]
                |> Array.ofList

let restructureParameterSet (input : parameterSearch) (oneDimensionalSurvival:float []) (oneDimensionalCloneSize: float [] []) =
    let rhoN = Array.length input.rhoRange
    let rN = Array.length input.rRange
    let lambdaN = Array.length input.lambdaRange
    let timePointsN = Array.length input.timePoints
    //Need to check the order of the matrix to avoid matrix transposition
    //printf "Length %A -> Expected %A\n" (Array.length oneDimensionalSurvival) (rhoN*rN*lambdaN*timePointsN)
    let (deltaRange,deltaN) = 
                    match input.deltaRange with
                    | Zero      -> ([0.],1)
                    | Range(r)  -> (List.map (fun i -> i*1.<probability^-1>) r,(List.length r))
    let probS = Array.init deltaN (fun delta -> Array.init lambdaN (fun lambda -> Array.init rhoN (fun rho -> Array.init rN (fun r -> Array.init timePointsN (fun t -> oneDimensionalSurvival.[t+(r*timePointsN)+(rho*timePointsN*rN)+(lambda*rhoN*timePointsN*rN)+(delta*lambdaN*rhoN*timePointsN*rN) ]  )))))
    let probN = Array.init deltaN (fun delta -> Array.init lambdaN (fun lambda -> Array.init rhoN (fun rho -> Array.init rN (fun r -> Array.init timePointsN (fun t -> oneDimensionalCloneSize.[t+(r*timePointsN)+(rho*timePointsN*rN)+(lambda*rhoN*timePointsN*rN)+(delta*lambdaN*rhoN*timePointsN*rN) ] )))))
    {   input with results = Some( 
                                    {
                                    scanResults.cloneSizeMatrix =   probN
                                    scanResults.survivalMatrix  =   probS
                                    scanResults.oneDimSizeMatrix =  oneDimensionalCloneSize
                                    scanResults.oneDimSurvMatrix =  oneDimensionalSurvival
                                    }    )
        }

let resultsMap f input =
    //Apply a function to every parameter in the results cloneSizeMatrix
    //For each parameter there is an array (time) of arrays (probability of seeing a number of clones, from 1 to maxN)
    //d * l * rho * r
    Array.map (fun d -> Array.map (fun l -> Array.map (fun rho -> Array.map (fun r -> f r) rho ) l ) d ) input

let resultsMapi f input =
    //Apply a function to every parameter in the results cloneSizeMatrix
    //For each parameter there is an array (time) of arrays (probability of seeing a number of clones, from 1 to maxN)
    //d * l * rho * r
    Array.mapi (fun id d -> Array.mapi (fun il l -> Array.mapi (fun irho rho -> Array.mapi (fun ir r -> f id il irho ir r) rho ) l ) d ) input

let resultsMap2 f input1 input2 =
    Array.map2 (fun d1 d2 -> Array.map2 (fun l1 l2 -> Array.map2 (fun rho1 rho2 -> Array.map2 (fun r1 r2 -> f r1 r2) rho1 rho2 ) l1 l2 ) d1 d2 ) input1 input2

let likelihoodTo1D input = 
    let dL = Array.length input
    let lL = Array.length input.[0]
    let rhoL = Array.length input.[0].[0]
    let rL = Array.length input.[0].[0].[0]
    Array.init (dL*lL*rhoL*rL) (fun i -> input.[(i/(rL*rhoL*lL))%dL].[(i/(rL*rhoL))%lL].[(i/rL)%rhoL].[i%rL])