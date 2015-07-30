module Whittaker

open MathNet.Numerics

let M k m (z:complex) = 
    let constant =  z**(m+0.5) * exp(-z/(complex 2. 0.)) 
    constant * (Kummer.M (m - k + 0.5) (1. + 2.*m) z) None

let W k m (z:complex) =
    let constant =  z**(m+0.5) * exp(-z/(complex 2. 0.)) 
    constant * (Kummer.U (m - k + 0.5) (1. + 2.*m) z) None
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*m) * (M k m z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - k - m) + MathNet.Numerics.SpecialFunctions.Gamma(2.*m) * (M k -m z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - k + m)

