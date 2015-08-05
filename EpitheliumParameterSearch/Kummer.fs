module Kummer

open MathNet.Numerics

let debug = true

//Todo: implement a gamma function which can accept complex arguements
//Note: 169 is the max input to this function before the value becomes infinite
let gamma = MathNet.Numerics.SpecialFunctions.Gamma

let rec complexGamma (z:System.Numerics.Complex) =
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
        let complexPi = complex System.Math.PI 0.
        if z.r < 0.5 then complexPi / ( (Trig.Sin (complexPi*z) ) * complexGamma ( (complex 1. 0.) - z) )
            else    let z = z - complex 1. 0.
                    //let x = complex 0.99999999999980993 0.
                    let x  = List.mapi (fun i pVal -> pVal/(z+complex (float(i)+1.) 0.)) p
                            |> List.fold (fun acc pVal -> acc+pVal ) (complex 0.99999999999980993 0.)  
                    let t = (complex (float(List.length p) - 0.5) 0.) + z
                    let result = sqrt( (complex 2. 0.)*complexPi) * t**(z+(complex 0.5 0.)) * exp(-t) * x
                    if debug then printf "complex gamma: z=%A result=%A\ngamma(%A+%Ai)\n" z result z.r z.i
                    result
    //gamma needs to be complex but lanczos algorithm as implemented is weak for small real z
    //Test and fallback to Numerics gamma if safe to do so
    //lanczos might work better if the function was implemented as a series computing to desired precision
    if z.i = 0. 
        then    let result = complex (gamma(z.r)) 0. 
                if debug then printf "Real gamma: z=%A result=%A\ngamma(%A)\n" z.r result z.r
                result
        else lanczos z


let factorial n =
    let rec core n acc =
        if n=0 then acc else core (n-1) (float(n)*acc)
    core n 1.

let hyperGeometric0F1 a z =
    //The confluent hypergeometric limit function
    //known weakness- large z
    let accuracy = pown 10. -15
    let rec core a z n delta result =
        let delta' = if n=0 then complex 1. 0. else delta * z / (complex ((float(n)+a-1.)*(float(n))) 0. )
        let result' = if n=0 then complex 1. 0. else result + delta'
        match( ((delta'/result').Magnitude < accuracy), n<500) with
        | (true,_) -> result'
        | (false,true) -> core a z (n+1) delta' result'
        | (false,false) -> failwith("The confluent hypergeometric limit function failed to converge")
    core a z 0 (complex 0. 0.) (complex 0. 0.)

let besselJ v x =
    (x/(complex 2. 0.))**(complex v 0.) / (complex (gamma(v+1.)) 0.) * (hyperGeometric0F1 (v+1.) (-x*x/(complex 4. 0.)) )

let M (a:complex) (b:complex) (z:complex) = 
    //Test the input for anything untoward- we cannot cope with NaN
    if System.Double.IsNaN(a.Magnitude) || System.Double.IsNaN(b.Magnitude) || System.Double.IsNaN(z.Magnitude) then failwith "Kummer M NaN input: a %A b %A z %A\n" a b z
    
    //Different ways to calculate the hypergeometric function, based on their strengths. 
    //Current known weaknesses
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
    let rec taylorExpansion accuracy a b z n numerator denominator result =
        //Taylor expansion of M. This has been reported as an efficent way of calculating M, but with known weaknesses
        //let numerator' = if n = 0 then complex 1. 0. else numerator * (complex (a + float(n-1)) 0.) * z
        //let denominator' = if n = 0 then complex 1. 0. else denominator * complex (float(n) * (b + float(n-1))) 0.
        let numerator' = if n = 0 then complex 1. 0. else numerator * (a + (complex (float(n-1)) 0.) ) * z
        let denominator' = if n = 0 then complex 1. 0. else denominator * (complex (float(n)) 0.)*(b + (complex (float(n-1)) 0.) )
        let delta = numerator' / denominator'
        let result' = result + delta
        if (delta/result).Magnitude < accuracy then ignore (if debug then printf "Taylor expansion in %A steps\na=%A b=%A c=%A result=%A\nhypergeom(%A+%Ai,%A+%Ai,%A+%Ai)\n" n a b z result' a.r a.i b.r b.i z.r z.i); result' else 
            if n+1 < 500 then taylorExpansion accuracy a b z (n+1) numerator' denominator' result' else printf "a %A b %A c %A r %A\n" a b z result' ; failwith("Kummer M function failed to converge (Taylor)")
    let rec singleFraction accuracy a b z n alpha beta gamma result result' =
        //Expressing the Taylor expansion (above) as a single fraction copes better when b < 1
        //let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * complex (float(n) * (b + float(n) - 1.)) 0.
        let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * (complex (float(n)) 0.) * (b + (complex (float(n)-1.) 0.) )  //complex (float(n) * (b + float(n) - 1.)) 0.
        let beta'  = if n = 0 then complex 1. 0. else beta * (a + complex (float(n) - 1.) 0.) * z
        let gamma' = if n = 0 then complex 1. 0. else gamma * (complex (float(n)) 0.) * (b + (complex (float(n) - 1.) 0.))

        let result'' = (alpha' + beta')/gamma'

        if (((result'' - result')/result').Magnitude < accuracy) && (((result' - result)/result).Magnitude < accuracy) then ignore (if debug then printf "Single fraction in %A steps\na=%A b=%A c=%A result=%A\nhypergeom(%A+%Ai,%A+%Ai,%A+%Ai)\n" n a b z result'' a.r a.i b.r b.i z.r z.i); result'' else
            if n+1 < 500 then singleFraction accuracy a b z (n+1) alpha' beta' gamma' result' result'' else  printf "a %A b %A c %A r %A\n" a b z result'; failwith("Kummer M function failed to converge (Fraction)")
    let rec buchholz accuracy a b z n d d' d'' result =
        //Initialise system
        //First three steps are skipped as we know the d coefficients prior to calculation
        let unchangingCoefficient = (complex (gamma(b) * 2. ** (b-1.)) 0.) * exp(z/(complex 2. 0.))
        let sqrt_z_x_4a_minus_2b = sqrt( z*(complex (4.*a-2.*b) 0.) )
        let result = if n = 0 then unchangingCoefficient * ( ( besselJ (b-1.)    ( sqrt_z_x_4a_minus_2b ) )/(sqrt_z_x_4a_minus_2b**( b-1.   )) +
                                                             ( besselJ (b-1.+2.) ( sqrt_z_x_4a_minus_2b ) )/(sqrt_z_x_4a_minus_2b**( b-1.+2. )) * (complex (b/2.) 0.) * z ** 2.  ) else result
        let n = if n = 0 then 3 else n
        let d''' = ( (complex (float(n) - 2. + b) 0.) * d' + (complex (2.*a - b) 0.) * d )/(complex (float(n)) 0.)
        let delta = unchangingCoefficient * d''' * (besselJ (b-1.+float(n)) (sqrt_z_x_4a_minus_2b ) ) / (sqrt_z_x_4a_minus_2b**( (b-1.+float(n)) ))
        let result' = result + delta
        match (((delta/result').Magnitude < accuracy) ,n<500) with 
        | (true,_) -> ignore (if debug then printf "Buchholz in %A steps\n" n); result'
        | (false,true) -> buchholz accuracy a b z (n+1) d' d'' d''' result'
        | (false,false) -> printf "a %A b %A c %A\n" a b z ; failwith("Kummer M function failed to converge (Buchholz)")

    
    //Decide which approach to use
    //Temporarily disabling Buchholz method until I've confirmed that I can reproduce the original matlab results
    //Need a test which copes with NaN sign(a)=sign(z.Real)
    //if debug then printf "a %A b %A z %A\n" a b z

    match (b.Magnitude>=1.,false) with 
    | (true, _) -> taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.)
    | (false, _) -> singleFraction accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 1. 0.) (complex 0. 0.) (complex 1. 0.)
    //| (_, false)    -> buchholz accuracy a b z 0 (complex 1. 0.) (complex 0. 0.) (complex 0. 0.) (complex 0. 0.)

let U a (b:complex) z = 
    //undefined for integer b, so we make small perturbations to integer
    let b = if b.r%1. = 0. then b + complex 0.00000001 0. else b
    //(M a b z)* (complex (gamma(1.-b)/gamma(1.+a-b)) 0. ) + (M (1.+a-b) (2.-b) z) * (complex (gamma(b-1.)/gamma(a)) 0.) * z**(1.-b)
    (M a b z)* ( complexGamma((complex 1. 0.)-b)/complexGamma((complex 1. 0.)+a-b) ) + (M ((complex 1. 0.)+a-b) ((complex 2. 0.)-b) z) * ( complexGamma(b-(complex 1. 0.))/complexGamma(a) ) * z**((complex 1. 0.)-b)
    
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*b) * (M a b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a - b) + MathNet.Numerics.SpecialFunctions.Gamma(2.*b) * (M a -b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a + b)