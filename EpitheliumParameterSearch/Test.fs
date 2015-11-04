module Test

open MathNet.Numerics

ignore (Kummer.M (complex 0.193813772152103 0.) (complex 0.999 0.) (complex 33.8845488347767 -10.6505935688991))


let gamma = 0.05263157895
let T =  0.1571428571<Types.cell>
let k_50 = (complex -0.0320515775716552 0.999486216200688)
let r = 0.01
let t = (T*1.<Types.cell^-1>) 
ignore (AnalyticalCloneSizeDistribution.F k_50 k_50 t r gamma)