module Types

open MathNet.Numerics

[<Measure>]
type cell

[<Measure>]
type week

[<Measure>]
type probability

type delta = Zero | Range of float<probability> list

type parameterSearch = {    rhoRange:       float array
                            rRange:         float<probability> array
                            lambdaRange:    float<cell/week> array
                            timePoints:     float<week> array
                            maxN:           int
                            deltaRange:     delta
                            cloneSizeMatrix:  float [] [] [] [] [] option // lambda x rho x r x time x n
                            survivalMatrix: float [] [] [] [] option //lambda x rho x r x time
                            oneDimSizeMatrix: float [] [] option
                            oneDimSurvMatrix: complex [] option
                            }

type parameterSet = {   time: float<week>;
                        rho: float; 
                        r: float<probability>;
                        lambda: float<cell/week>
                        } with 
                        member this.migration = this.rho/(1.-this.rho) * this.lambda
                        member this.correctPathologicalPoint = { this with rho = this.rho + 0.0001 ; r = this.r+0.00001<probability> }
                        override this.ToString() = sprintf "Rho: %A r:%A Time:%A Lambda:%A" this.rho this.r this.time this.lambda

let testSystem = {time=1.<week>; rho=0.85; r=0.15<probability>; lambda=2.<cell/week>}
