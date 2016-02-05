﻿module GeneratePlot

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

let maxLikelihood (matrix: float [] [] [] []) = 
    let mutable max = -infinity
    let mutable coordinate = (0,0,0,0)
    let a = Array.length matrix
    let b = Array.length matrix.[0]
    let c = Array.length matrix.[0].[0]
    let d = Array.length matrix.[0].[0].[0]

    for i=0 to (a-1) do
        for j=0 to (b-1) do
            for k=0 to (c-1) do
                for l=0 to (d-1) do
                    if matrix.[i].[j].[k].[l] > max then 
                        max <- matrix.[i].[j].[k].[l]
                        coordinate <- (i,j,k,l)

    (max,coordinate)

let fourDArrayMap (matrix: float [] [] [] []) f = 
    Array.map (fun i2 -> Array.map (fun i1 -> Array.map (fun i0 -> Array.map f i0) i1) i2) matrix

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
    file.Close

//let projectDeltaXRho likelihoodMatrix = 
    


