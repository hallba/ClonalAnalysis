module Gamma

open MathNet.Numerics

let debug = true

let eulerConstant = 0.5772156649015329
let eulerComplex  = complex eulerConstant 0.

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
    if debug then printf "Gamma %A = %A\n" z result
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
 
let diGammaInt x = 
    let rec core acc n =
        let n' = n - 1
        if n' = 0 then acc else core (acc+1./float(n')) n'
    if x <= 0 then infinity else (core (0.-eulerConstant) x)

//Coefficients for C.Lanczos expansion of DiGamma function dg_coef[k] = - (2k+1) * lg_coef[k]
let dg_coeff = [|   -0.83333333333333333e-1; 
                    0.83333333333333333e-2;
                    -0.39682539682539683e-2;
                    0.41666666666666667e-2;
                    -0.75757575757575758e-2;
                    0.21092796092796093e-1;
                    -0.83333333333333333e-1;
                    0.4432598039215686;
                    -0.3053954330270122e+1; 
                    0.125318899521531e+2;
                    |]
//-1/12, 3/360,-5/1260, 7/1680,-9/1188, 11*691/360360,-13/156, 15*3617/122400, ? , ?

let rec diGammaFloat x =
    if x < 0. then (diGammaFloat (1.-x) - System.Math.PI/tan(System.Math.PI*x)) else
        let xa = abs(x)
        if x%1. = 0. then diGammaInt (int(x)) 
        else if (xa+0.5)%1. = 0. && x > 0.5 then 
            let n = int (xa - 0.5) 
            List.init n (fun k -> float(k+k+1))
            |> List.fold (fun acc n -> acc + 1./n ) 0.
            |> fun dgam -> 2.*dgam - 2.*log(2.) - eulerConstant
        else 
            let (dgam, xa) =    if xa >= 10. then (0.,xa) 
                                else List.init (10 - int(xa)) (fun i -> i)
                                     |> List.fold (fun acc k -> acc - 1./(float(k)+xa) ) 0.
                                     |> fun dgam -> (dgam, (xa+ 10. - floor(xa)) )
            let dgam' = dgam + log(xa) - 0.5/xa
            let overx2 = 1./(xa*xa)
            List.init 10 (fun i -> ( i, (pown overx2 (i+1)) ) )
            |> List.fold (fun acc (k,ov) -> acc + dg_coeff.[k]*ov ) dgam'

let diGammaComplex (z:complex) =
    if z.i = 0. && z.r % 1. <> 0. then (complex (diGammaFloat z.r) 0.)  //x is actually real and a float
    else if z.i = 0. then (complex (diGammaInt (int(z.r))) 0.)          //x is actually real and an int
    else    let (zp,zra) = if z.r < 0. then (-z,-z.r)  else (z,z.r)
            let n = 0
            let (zm,n) = if zra<8. then ((zp+(complex (float(8-int(z.r))) 0.)),(8-int(z.r))) else (zp,n)
            let overz = (complex 1. 0.) / zm
            let overz2 = overz*overz
            let dgam =          log(zm) - overz/(complex 2. 0.)
            let dgam' =         List.init 10 (fun k -> ( k, (overz2**2.) ) )
                                |> List.fold (fun acc (k,ok) -> acc + (complex dg_coeff.[k] 0.)*ok) dgam
            let dgam''=         List.init 10 (fun k -> zp + (complex (float(k)) 0.)*(complex 1. 0.) )
                                |> List.fold (fun acc k -> acc - (complex 1. 0.)/k ) dgam'
            dgam''