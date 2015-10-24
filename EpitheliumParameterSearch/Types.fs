module Types

open MathNet.Numerics

[<Measure>]
type cell

[<Measure>]
type week

[<Measure>]
type probability

type parameterSearch = {    rhoRange:       float array
                            rRange:         float<probability> array
                            lambdaRange:    float<cell/week> array
                            timePoints:     float<week> array
                            maxN:           int

                            cloneSizeMatrix:  float [] [] [] [] [] option // lambda x rho x r x time x n
                            survivalMatrix: float [] [] [] [] option //lambda x rho x r x time
                            oneDimSizeMatrix: float [] [] option
                            oneDimSurvMatrix: complex [] option
                            }