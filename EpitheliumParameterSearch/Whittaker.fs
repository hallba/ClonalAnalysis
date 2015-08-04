﻿module Whittaker

open MathNet.Numerics

let M (k:complex) m (z:complex) = 
    //let constant =  z**(complex (m+0.5) 0. ) * exp(-z/(complex 2. 0.)) 
    let constant =  z**( m + (complex 0.5 0.) ) * exp(-z/(complex 2. 0.)) 
    if System.Double.IsNaN(k.Magnitude) || System.Double.IsNaN(m.Magnitude) || System.Double.IsNaN(z.Magnitude)  then failwith "Whittaker M NaN input: a %A b %A z %A\n" m k z
    constant * (Kummer.M (m - k + (complex 0.5 0.) ) ( (complex 1. 0.) + (complex 2. 0.)*m) z) 

let W (k:complex) m (z:complex) =
    let constant =  z**(m + complex 0.5 0. ) * exp(-z/(complex 2. 0.)) 
    if System.Double.IsNaN(k.Magnitude) || System.Double.IsNaN(m.Magnitude) || System.Double.IsNaN(z.Magnitude) then failwith "Whittaker W NaN input: a %A b %A z %A\n" m k z
    constant * (Kummer.U (m - k + (complex 0.5 0.) ) ( (complex 1. 0.) + (complex 2. 0.)*m) z) 
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*m) * (M k m z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - k - m) + MathNet.Numerics.SpecialFunctions.Gamma(2.*m) * (M k -m z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - k + m)