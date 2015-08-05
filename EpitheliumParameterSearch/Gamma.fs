module Gamma

open MathNet.Numerics

let debug = true

//Todo: implement a gamma function which can accept complex arguements
//Note: 169 is the max input to this function before the value becomes infinite
let realGamma = MathNet.Numerics.SpecialFunctions.Gamma

let logGamma (z:complex) =
    //Stirlings formula
    let coefficients = [    complex (1./12.) 0.;
                            complex (-1./360.) 0.;
                            complex (1./1260.) 0.;
                            complex (-1./1680.) 0.;
                            complex (1./1188.) 0.;
                            complex (-691./360360.) 0.;
                            complex (1./156.) 0.;
                            complex (-3617./122400.) 0.;
                            complex (43867./244188.) 0.;
                            complex -1.39243221690590 0.  
                            ]
    let result = List.mapi (fun i c -> c / (z**(complex (float(i*2+1)) 0. )) ) coefficients
                    |> List.fold (fun acc item -> acc + item) ( (z - (complex 0.5 0.))*log(z) - z + (complex (0.5*log(2.*System.Math.PI)) 0.) )
    result

let stirling z =
    exp (logGamma z)

let complexPi = complex System.Math.PI 0.

let reflection nextFunction z =
    let oneMinusZ = complex 1. 0. - z
    complexPi / ( (Trig.Sin (complexPi*z) ) * nextFunction ( oneMinusZ ) )       

let lanczos (z:complex) =
    //Lanczos approximation, based on wikipedia implementation
    //http://creativecommons.org/licenses/by-sa/3.0/
    let p = [   complex 676.5203681218851 0.;
                complex -1259.1392167224028 0.;
                complex 771.32342877765313 0.;
                complex -176.61502916214059 0.;
                complex 12.507343278686905 0.;
                complex -0.13857109526572012 0.;
                complex 9.9843695780195716e-6 0.;
                complex 1.5056327351493116e-7 0.  ]
    let z' = z - complex 1. 0.
                //let x = complex 0.99999999999980993 0.
    let x  = List.mapi (fun i pVal -> pVal/(z'+complex (float(i)+1.) 0.)) p
                        |> List.fold (fun acc pVal -> acc+pVal ) (complex 0.99999999999980993 0.)  
    let t = (complex (float(List.length p) - 0.5) 0.) + z'
    let result = sqrt( (complex 2. 0.)*complexPi) * t**(z'+(complex 0.5 0.)) * exp(-t) * x
    result

let rec complexGamma (z:System.Numerics.Complex) =
    //gamma needs to be complex but Stirlings/Lanczos are less precise approximations
    //Better to keep to one function or mix/match? Stirling gives NaN for too many real z but lanczos is too imprecise
    if z.i = 0. 
        then    let result = complex (realGamma(z.r)) 0. 
                if debug then printf "Real gamma: z=%A result=%A\ngamma(%A)\n" z.r result z.r
                result
        elif z.r < 0.5 then reflection lanczos z
        else
                let result = lanczos z
                if debug then printf "Complex gamma: z=%A result=%A\ncgama(%A,%A,1)\n" z result z.r z.i
                result