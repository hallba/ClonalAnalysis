module IO

let storeSearch (filename: string) (search: Types.parameterSearch) = 
        use file = new System.IO.FileStream(filename, System.IO.FileMode.Create)
        let bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter()
        use mem = new System.IO.MemoryStream()
        bf.Serialize(mem,search)
        let byteArray = mem.ToArray()
        file.Write(byteArray,0,byteArray.Length)
        file.Close()
        ()

let retrieveSearch (filename: string) =
    use fileStream = new System.IO.FileStream(filename, System.IO.FileMode.Open)
    let bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter()
    let result = bf.Deserialize(fileStream)
    fileStream.Close()
    (result :?> Types.parameterSearch)

//Get data from a matlab generated matrix
let importMatlab (filename: string) = 

    let m = new csmatio.io.MatFileReader(filename)
    
    let read1d name =
        let r = m.Content.[name] :?> csmatio.types.MLDouble
        let s = r.Dimensions
        let i = r.GetArray()
        i.[0]

    //readNd - Read a multidimensional matrix from an open matlab file and return an array of arrays (of arrays of arrays...)
    let readNd name (dimensions:int list) = 
        let r = m.Content.[name] :?> csmatio.types.MLDouble

        //First- read all of the values into a 1D array
        let rec totalElements d acc = 
            match d with 
            | [] -> acc
            | topD::rest -> totalElements rest (acc*topD)
        
        Array.init (totalElements dimensions 1) (fun i -> r.Get(i))

    let readCells name = 
        let r = m.Content.[name] :?> csmatio.types.MLCell
        Array.init r.Size (fun i -> (r.Cells.[i] :?> csmatio.types.MLDouble).GetArray())

    let rRange = read1d "rRange" |> Array.map (fun i -> i*1.<Types.probability>)
    let lambdaRange = read1d "lambdaRange" |> Array.map (fun i -> i*1.<Types.cell/Types.week>)
    let rhoRange = read1d "rhoRange" |> Array.map (fun i -> 1. - i) //Rho correction due to gamma inversion in matlab code
    let t = read1d "t" |> Array.map (fun i -> i*1.<Types.week>)
    let timePoints = read1d "timePoints"
    let nRange = readCells "nRange"
    let maxN = int (read1d "maxN").[0]

    let input = {   Types.rhoRange=rhoRange
                    Types.rRange=rRange
                    Types.lambdaRange=lambdaRange
                    Types.deltaRange=Types.Zero
                    Types.timePoints=t
                    Types.maxN=maxN
                    Types.results=None
                    Types.excludeOnes=true
                    }

    let cloneSizeP = readNd "PScanPP"  [(Array.length lambdaRange);(Array.length rhoRange);(Array.length rRange);(Array.length t)] //This is a problem- n is of variable length
    let survivalP  = readNd "PSurvScanPP" [(Array.length lambdaRange);(Array.length rhoRange);(Array.length rRange);(Array.length t)]

    Types.restructureParameterSet input survivalP cloneSizeP

    //could also read maxN and nBadvalues

//    let scanMatD = ([lambdaRange;rhoRange;rRange;t;nRange] |> List.map (fun i -> Array.length i))
//    let cloneSizeP = readNd "PScanPP" scanMatD
//                     |> fun matrix ->
//                        Array.init scanMatD.[0] 
//                            (fun i -> (Array.init scanMatD.[1] 
//                                (fun j -> (Array.init scanMatD.[2] 
//                                    (fun k -> ( Array.init scanMatD.[3]
//                                        (fun l -> ( Array.init scanMatD.[4] 
//                                            (fun m -> matrix.[i+(j*scanMatD.[0])+(k*scanMatD.[0]*scanMatD.[1])+(l*scanMatD.[0]*scanMatD.[1]*scanMatD.[2])+(m*scanMatD.[0]*scanMatD.[1]*scanMatD.[2]*scanMatD.[3])] 
//                                                ) ) ) ) ) ) ) ) )
//    let survivalP  = readNd "PSurvScanPP" ([lambdaRange;rhoRange;rRange;t] |> List.map (fun i -> Array.length i))
