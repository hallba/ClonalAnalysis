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




