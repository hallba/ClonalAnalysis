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

//Note- Stirling is poor for small z
let stirling z =
    exp (logGamma z)

let complexPi = complex System.Math.PI 0.

let reflection nextFunction z =
    let oneMinusZ = complex 1. 0. - z
    complexPi / ( (Trig.Sin (complexPi*z) ) * nextFunction ( oneMinusZ ) )       

let lanczos (z:complex) =
    //Lanczos approximation, based on wikipedia 
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

let lanczosGodfrey (z:complex) =
//The "best" set of found coefficients without swapping for a polynomial expansion, ala Boost
//From Godfrey http://my.fit.edu/~gabdo/gamma.txt
//                          (z + 0.5)
//             (z + g + 0.5)
//Gamma(z+1) = (-----------)        * S
//             (    e      )
//
//
//with
//
//             a(1)    a(2)    a(3)    a(4)    a(5)    a(6)
//S ~= [a(0) + ----- + ----- + ----- + ----- + ----- + ----- + ... ]
//             z + 1   z + 2   z + 3   z + 4   z + 5   z + 6
    let c = [  0.99999999999999709182 ;
          57.156235665862923517 ;
         -59.597960355475491248 ;
          14.136097974741747174 ;
          -0.49191381609762019978 ;
           0.33994649984811888699e-4 ;
           0.46523628927048575665e-4 ;
          -0.98374475304879564677e-4 ;
           0.15808870322491248884e-3 ;
          -0.21026444172410488319e-3 ;
           0.21743961811521264320e-3 ;
          -0.16431810653676389022e-3 ;
           0.84418223983852743293e-4 ;
          -0.26190838401581408670e-4 ;
           0.36899182659531622704e-5 ] 
    let g = complex (607./128.) 0.
    let W = exp(g)/sqrt(complex (2.*System.Math.PI) 0.)
    let c0_5 = complex 0.5 0.
    let z' = z - complex 1. 0.
    let ce = complex System.Math.E 0.
    //let a = [0.01689 ; 1.2866 ; -1.461 ; 0.4055 ; -0.02080 ; 2.0413E-05 ; -9.1123E-08]
    //let g = complex 5. 0.
    let result = List.map (fun f -> complex f 0.) c
                 |> List.map (fun c -> c/W ) //Skip this step if using a above
                 |> List.mapi (fun i c -> if i = 0 then c else c/(z'+(complex (float(i)) 0. ) ) )
                 |> List.fold (fun acc item -> acc + item ) (complex 0. 0.)
                 |> fun f ->  ( ( (z' + g + c0_5 )/ce ) ** (z'+c0_5) ) * f
    result
 
let logLanczosGodfrey (z:complex) =
//The "best" set of found coefficients without swapping for a polynomial expansion, ala Boost
//From Godfrey http://my.fit.edu/~gabdo/gamma.txt
//                          (z + 0.5)
//             (z + g + 0.5)
//Gamma(z+1) = (-----------)        * S
//             (    e      )
//
//
//with
//
//             a(1)    a(2)    a(3)    a(4)    a(5)    a(6)
//S ~= [a(0) + ----- + ----- + ----- + ----- + ----- + ----- + ... ]
//             z + 1   z + 2   z + 3   z + 4   z + 5   z + 6
    let c = [  0.99999999999999709182 ;
          57.156235665862923517 ;
         -59.597960355475491248 ;
          14.136097974741747174 ;
          -0.49191381609762019978 ;
           0.33994649984811888699e-4 ;
           0.46523628927048575665e-4 ;
          -0.98374475304879564677e-4 ;
           0.15808870322491248884e-3 ;
          -0.21026444172410488319e-3 ;
           0.21743961811521264320e-3 ;
          -0.16431810653676389022e-3 ;
           0.84418223983852743293e-4 ;
          -0.26190838401581408670e-4 ;
           0.36899182659531622704e-5 ] 
    let g = complex (607./128.) 0.
    let W = exp(g)/sqrt(complex (2.*System.Math.PI) 0.)
    let c0_5 = complex 0.5 0.
    let z' = z - complex 1. 0.
    let ce = complex System.Math.E 0.
    //let a = [0.01689 ; 1.2866 ; -1.461 ; 0.4055 ; -0.02080 ; 2.0413E-05 ; -9.1123E-08]
    //let g = complex 5. 0.
    let result = List.map (fun f -> complex f 0.) c
                 |> List.map (fun c -> c/W ) //Skip this step if using a above
                 |> List.mapi (fun i c -> if i = 0 then c else c/(z'+(complex (float(i)) 0. ) ) )
                 |> List.fold (fun acc item -> acc + item ) (complex 0. 0.)
                 |> fun f ->  (log f) - z' - c0_5 + (z'+c0_5)*log(z'+g+c0_5) 
    result
    

let rec complexGamma algorithm (z:System.Numerics.Complex) =
    //gamma needs to be complex but Stirlings/Lanczos are less precise approximations
    //Better to keep to one function or mix/match? Stirling gives NaN for too many real z but lanczos is too imprecise
    if z.i = 0. 
        then    let result = complex (realGamma(z.r)) 0. 
                if debug then printf "Real gamma: z=%A result=%A\ngamma(%A)\n" z.r result z.r
                result
        elif z.r < 0.5 then reflection algorithm z
        else
                let result = algorithm z
                if debug then printf "Complex gamma: z=%A result=%A\ncgama(%A,%A,1)\n" z result z.r z.i
                result