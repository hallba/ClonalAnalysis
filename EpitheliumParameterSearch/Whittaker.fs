module Whittaker

open MathNet.Numerics

let M k m z = 
    let constant =  z**(complex (m+0.5) 0. ) * exp(-z/(complex 2. 0.)) 
    if System.Double.IsNaN(k) || System.Double.IsNaN(m) || System.Double.IsNaN(z.r) || System.Double.IsNaN(z.i) then failwith "Whittaker M NaN input: a %A b %A z %A\n" m k z
    constant * (Kummer.M (m - k + 0.5) (1. + 2.*m) z) 

let W k m z =
    let constant =  z**(complex (m+0.5) 0. ) * exp(-z/(complex 2. 0.)) 
    if System.Double.IsNaN(k) || System.Double.IsNaN(m) || System.Double.IsNaN(z.r) || System.Double.IsNaN(z.i) then failwith "Whittaker W NaN input: a %A b %A z %A\n" m k z
    constant * (Kummer.U (m - k + 0.5) (1. + 2.*m) z) 
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*m) * (M k m z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - k - m) + MathNet.Numerics.SpecialFunctions.Gamma(2.*m) * (M k -m z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - k + m)

