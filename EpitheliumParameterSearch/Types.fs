module Types

open MathNet.Numerics

[<Measure>]
type cell

[<Measure>]
type week

[<Measure>]
type probability

type delta = Zero | Range of float<probability> list

type scanResults       = {  cloneSizeMatrix:  float [] [] [] [] []  // lambda x rho x r x time x n
                            survivalMatrix: float [] [] [] []  //lambda x rho x r x time
                            oneDimSizeMatrix: float [] [] 
                            oneDimSurvMatrix: complex [] 
                            }

type parameterSearch = {    rhoRange:       float array
                            rRange:         float<probability> array
                            lambdaRange:    float<cell/week> array
                            timePoints:     float<week> array
                            maxN:           int
                            deltaRange:     delta
                            results:        scanResults option
                            }

type parameterSet = {   time: float<week>;
                        rho: float; 
                        r: float<probability>;
                        delta: float<probability>
                        lambda: float<cell/week>
                        } with 
                        member this.migration = this.rho/(1.-this.rho) * this.lambda
                        member this.correctPathologicalPoint = { this with rho = this.rho + 0.0001 ; r = this.r+0.00001<probability> }
                        override this.ToString() = sprintf "Rho: %A r:%A Time:%A Lambda:%A" this.rho this.r this.time this.lambda

let testSystem = {time=1.<week>; rho=0.85; r=0.15<probability>; delta=0.<probability>; lambda=2.<cell/week>}

let createParameterSet (input : parameterSearch) = 
    [ for lambda in input.lambdaRange do for rho in input.rhoRange do for r in input.rRange do for t in input.timePoints do yield {testSystem with lambda=lambda;time=t;r=r;rho=rho} ]
    |> Array.ofList