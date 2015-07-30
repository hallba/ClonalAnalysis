module Kummer

open MathNet.Numerics

let debug = true

let gamma = MathNet.Numerics.SpecialFunctions.Gamma

let factorial n =
    let rec core n acc =
        if n=0 then acc else core (n-1) (float(n)*acc)
    core n 1.

let besselJ v x =
    let accuracy = pown 10. -15
    let rec core v x n result = 
        //let delta = (pown -1. n) * ((x/2.)**(2.*(float(n) + v))) / ( (factorial n) * gamma(v+float(n)+1.) )
        let delta = (complex (pown -1. n) 0.) / (complex ( (factorial n) * gamma(v+float(n)+1.) ) 0.) * ( x/(complex 2. 0.) ) ** (2.*(float(n) + v)) 
        let result' = result + delta
        if (delta/result).Magnitude < accuracy then result else 
            if n < 500 then core v x (n+1) result' else printf "Last result = %A\n" result'; failwith("Bessel J function failed to converge")
    core v x 0 (complex 0. 0.)
    

let M a b (z:complex) maxN = 
    //Different ways to calculate the hypergeometric function, based on their strengths. 
    //Current known weaknesses
    //-Anything with complex a or b
    //--because of the U function below- which requires a gamma function which handles complex numbers
    //-z > 100
    //-a or b > 50
    //-sign(Re(a)) = negative sign(Re(z))
    //
    //Things to implement in future
    //-Recurrance relation ((z> 100 or (a>50 or b>50))
    //-Asymptotic (z>100 and a<50 and b<50)
    //-Bucholz polynomial (sign(a) = negative sign(z), a<50, b<50, z<100
    //
    //See Masters Thesis of John Pearson for more details https://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf
    let accuracy = pown 10. -15
    let maxAttempts = match maxN with
                        | Some(n) -> n
                        | None -> 500
    let rec taylorExpansion accuracy a b z n numerator denominator result =
        //Taylor expansion of M. This has been reported as an efficent way of calculating M, but with known weaknesses
        let numerator' = if n = 0 then complex 1. 0. else numerator * (complex (a + float(n-1)) 0.) * z
        let denominator' = if n = 0 then complex 1. 0. else denominator * complex (float(n) * (b + float(n-1))) 0.
        let delta = numerator' / denominator'
        let result' = result + delta
        //printf "N %A D %A R %A A %A\n" numerator' denominator' result' (result'-result)
        if (delta/result).Magnitude < accuracy then result' else 
            if n+1 < 500 then taylorExpansion accuracy a b z (n+1) numerator' denominator' result' else printf "a %A b %A c %A\n" a b z ; failwith("Kummer M function failed to converge (Taylor)")
    let rec singleFraction accuracy a b z n alpha beta gamma result result' =
        //Expressing the Taylor expansion (above) as a single fraction copes better when b < 1
        let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * complex (float(n) * (b + float(n) - 1.)) 0.
        let beta'  = if n = 0 then complex 1. 0. else beta * (complex (a + float(n) - 1.) 0.) * z
        let gamma' = if n = 0 then complex 1. 0. else gamma * (complex (float(n) * (b + float(n) - 1.)) 0.)

        let result'' = (alpha' + beta')/gamma'

        if (((result'' - result')/result').Magnitude < accuracy) && (((result' - result)/result).Magnitude < accuracy) then result'' else
            if n+1 < 500 then singleFraction accuracy a b z (n+1) alpha' beta' gamma' result' result'' else  printf "a %A b %A c %A\n" a b z ; failwith("Kummer M function failed to converge (Fraction)")
    let rec buchholz accuracy a b z n d d' d'' constant result =
        let constant = if n=0 then (complex (gamma(b)) 0.) * exp(z/(complex 2. 0.)) else constant
        let d''' = match n with
                    | 0 -> complex 1. 0.
                    | 1 -> complex 0. 0.
                    | 2 -> complex (b/2.) 0.
                    | _ -> (complex (float(n) - 2. + b) 0.) * d' + (complex (2.*a - b) 0.) * d
        //let delta = constant * d''' * z ** (complex (float(n)) 0.) * (complex ( besselJ ( b-1.+float(n)) ) 0.) * sqrt(z*(complex (2.*b-4.*a) 0.)) / ( (z * (complex (2.*b-4.*a) 0.) )**( 0.5 * (b-1.+float(n)) ) ) 
        
        let delta = constant * d''' * z ** (float(n)) * (besselJ (b-1.+float(n)) ( sqrt(z*(complex (2.*b - 4.*a) 0.) ) ) ) / (z*(complex (2.*b - 4.*a) 0.))**(0.5*(b-1.+float(n)))
        let result' = result + delta
        
        ignore ((delta/result').Magnitude < accuracy)

        //This function is incomplete.
        failwith("The buchholz method has not been completely implemented")
    
    //Decide which approach to use
    //sign(a)=sign(z.Real)
    if debug then printf "a %A b %A z %A\n" a b z

    match (b>=1.,true) with 
    | (true, true) -> taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.)
    | (false, true) -> singleFraction accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 1. 0.) (complex 0. 0.) (complex 1. 0.)
    | (_, false)    -> buchholz accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 0. 0.) (complex 0. 0.) (complex 0. 0.)
//    if b > 1. then taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.) else
//        singleFraction accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 1. 0.) (complex 0. 0.) (complex 1. 0.)

let U a b z maxN = 
    //undefined for integer b, so we make small perturbations to integer
    let b = if b%1. = 0. then b + 0.00000001 else b
    (M a b z maxN)* (complex (gamma(1.-b)/gamma(1.+a-b)) 0. ) + (M (1.+a-b) (2.-b) z maxN) * (complex (gamma(b-1.)/gamma(a)) 0.) * z**(1.-b)
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*b) * (M a b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a - b) + MathNet.Numerics.SpecialFunctions.Gamma(2.*b) * (M a -b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a + b)