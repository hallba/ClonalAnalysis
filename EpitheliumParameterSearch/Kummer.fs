module Kummer

open MathNet.Numerics

let gamma = MathNet.Numerics.SpecialFunctions.Gamma

let M a b z maxN = 
    //Different ways to calculate the hypergeometric function, based on their strengths. 
    //Current known weaknesses
    //-Anything with complex a or b
    //--only because of the U function below- which requires a gamma function which handles complex numbers
    //-z > 100
    //-a or b > 50
    //-sign(a) = negative sign(z)
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
            if n+1 < 500 then taylorExpansion accuracy a b z (n+1) numerator' denominator' result' else failwith("Kummer M function failed to converge")
    let rec singleFraction accuracy a b z n alpha beta gamma result result' =
        //Expressing the Taylor expansion (above) as a single fraction copes better when b < 1
        let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * complex (float(n) * (b + float(n) - 1.)) 0.
        let beta'  = if n = 0 then complex 1. 0. else beta * (complex (a + float(n) - 1.) 0.) * z
        let gamma' = if n = 0 then complex 1. 0. else gamma * (complex (float(n) * (b + float(n) - 1.)) 0.)

        let result'' = (alpha' + beta')/gamma'

        if (((result'' - result')/result').Magnitude < accuracy) && (((result' - result)/result).Magnitude < accuracy) then result'' else
            if n+1 < 500 then singleFraction accuracy a b z (n+1) alpha' beta' gamma' result' result'' else failwith("Kummer M function failed to converge")

    if b > 1. then taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.) else
        singleFraction accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 1. 0.) (complex 0. 0.) (complex 1. 0.)

let U a b z maxN = 
    //undefined for integer b, so we make small perturbations to integer
    let b = if b%1. = 0. then b + 0.00000001 else b
    (M a b z maxN) * (complex (gamma(1.-b) / gamma(1.+a-b)) 0. ) + (M (1.+a-b) (2.-b) z maxN) * (complex (gamma(b-1.)/gamma(a)) 0.) * z**(1.-b)
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*b) * (M a b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a - b) + MathNet.Numerics.SpecialFunctions.Gamma(2.*b) * (M a -b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a + b)