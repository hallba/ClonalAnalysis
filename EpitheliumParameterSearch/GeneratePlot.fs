module GeneratePlot

open System
open System.Collections.Generic
open System.IO

let twoDArrayToString matrix =
    let xL = Array.length matrix
    let yL = Array.length matrix.[0]
    let mutable s = ""
    for i=0 to (xL - 1) do
        for j=0 to (yL - 1) do
            s <- s + (sprintf " %f" matrix.[i].[j])
        s <- s + "\n"
    s

let compLikelihood (matrix: float [] [] [] []) f = 
    let mutable max = matrix.[0].[0].[0].[0]
    let mutable coordinate = (0,0,0,0)
    let a = Array.length matrix
    let b = Array.length matrix.[0]
    let c = Array.length matrix.[0].[0]
    let d = Array.length matrix.[0].[0].[0]

    for i=0 to (a-1) do
        for j=0 to (b-1) do
            for k=0 to (c-1) do
                for l=0 to (d-1) do
                    if (f matrix.[i].[j].[k].[l] max) then 
                        max <- matrix.[i].[j].[k].[l]
                        coordinate <- (i,j,k,l)

    (max,coordinate)

let maxLikelihood matrix =
    compLikelihood matrix (fun a b -> a > b)

let minLikelihood matrix = 
    compLikelihood matrix (fun a b -> a < b)
let fourDArrayMap (matrix: float [] [] [] []) f = 
    Array.map (fun i2 -> Array.map (fun i1 -> Array.map (fun i0 -> Array.map f i0) i1) i2) matrix

type heatmapText = {    scaleName:  string
                        title:      string
                        colorMap:   string option
                        }

let twoDArrayToGnuplot matrix filename text =
    use file = new StreamWriter(filename, false)
    file.WriteLine(sprintf "set title %A" text.title)
    file.WriteLine("unset key")
    file.WriteLine("set tic scale 0")
    ignore( 
        match text.colorMap with
        | None -> file.WriteLine("set palette defined (0 0 0 0, 1 0.5 0 0, 2 1 0 0, 3 1 0.5 0, 4 1 1 0, 5 1 1 1)") 
        | Some(p) -> file.WriteLine(sprintf "set palette defined %A" p)
        )
    file.WriteLine("$map1 << EOD")
    file.WriteLine(twoDArrayToString matrix)
    file.WriteLine("EOD")
    file.Close()
    ()

let twoDArrayToXfarbe matrix filename (title: string option) = 
    use file = new StreamWriter(filename, false)
    let xL = Array.length matrix
    let yL = Array.length matrix.[0]
    ignore (match title with
            | None -> ()
            | Some(s) -> file.WriteLine(s)
            )
    file.WriteLine(sprintf "%A %A" xL yL)
    file.WriteLine(twoDArrayToString matrix)
    file.Close()
    ()

let projectDeltaXRho (likelihoodMatrix:float [] [] [] []) = 
    let dL = Array.length likelihoodMatrix
    let rhoL = Array.length likelihoodMatrix.[0].[0]

    let lambdaL = Array.length likelihoodMatrix.[0]
    let rL = Array.length likelihoodMatrix.[0].[0].[0]

    let result = Array.init dL ( fun i -> Array.init rhoL (fun j -> 0.) )

    for i=0 to dL-1 do
        for j=0 to rhoL-1 do
            
            for k=0 to lambdaL-1 do
                for l=0 to rL-1 do
                    result.[i].[j] <- result.[i].[j] + likelihoodMatrix.[i].[k].[j].[l]

    result