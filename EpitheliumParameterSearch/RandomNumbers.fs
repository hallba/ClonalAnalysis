module RandomNumbers

let rec gaussianMP : System.Random -> float -> float -> float =
    //Random number generator based on Marsaglia Polar method
    //This calculates two numbers at a time so one is stored in a ref
    let stored = ref None
    (fun rng mean sd->
        //printf "mean %f sd %f\n" mean sd
        match !stored with 
        | Some(v) ->    //printf "Using stored value\n"
                        stored := None
                        v*sd+mean
        | None ->       //printf "Calculating two new values\n"
                        let u = rng.NextDouble()*2.-1.
                        let v = rng.NextDouble()*2.-1.
                        let w = u*u + v*v
                        //printf "u %f v %f w %f\n" u v w
                        if w < 1. then 
                            let mul = sqrt(-2.*log(w)/w)
                            //printf "mul %f" mul
                            stored := Some(v*mul)
                            u*mul*sd+mean
                        else
                            //Error in calculation
                            gaussianMP rng mean sd
    )   