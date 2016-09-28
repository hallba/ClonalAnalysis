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

let rec totalElements d acc = 
    match d with 
    | [] -> acc
    | topD::rest -> totalElements rest (acc*topD)

//Get data from a matlab generated matrix
let importMatlab (filename: string) = 

    let m = new csmatio.io.MatFileReader(filename)
    
    let read1d name =
        let r = m.Content.[name] :?> csmatio.types.MLDouble
        let s = r.Dimensions
        let i = r.GetArray()
        i.[0]

    //readNd - Read a multidimensional matrix from an open matlab file and return an array of arrays (of arrays of arrays...)
//    let readNd name (dimensions:int list) = 
//        //Not clear how to write arbitrary dimensions
//        let r = m.Content.[name] :?> csmatio.types.MLDouble
//        //This is ordered in the opposite way from previous 1D arrays in this code. 
//        let data = Array.init (totalElements dimensions 1) (fun i -> r.Get(i))
//        let rec core dim acc f =
//            match dim with
//            | [] -> acc
//            | topD::tail -> let acc' = Array.init topD 
//                            core tail (acc acc')
//        data

    let read4d name (dimensions:int list) = 
        if List.length dimensions <> 4 then failwith "read4d requires 4 dimensions"
        let r = m.Content.[name] :?> csmatio.types.MLDouble
        //This is ordered in the opposite way from previous 1D arrays in this code. 
        let data = Array.init (totalElements dimensions 1) (fun i -> r.Get(i))
        
        let dataResult = Array.init dimensions.[0] (
                            fun i -> (Array.init dimensions.[1] (
                                fun j -> (Array.init dimensions.[2] (
                                    fun k -> (Array.init dimensions.[3] (
                                        fun l -> data.[i+dimensions.[0]*j+dimensions.[1]*dimensions.[0]*k+dimensions.[2]*dimensions.[1]*dimensions.[0]*l] ) ) ) ) ) ) )

        dataResult

    let read5d name (dimensions:int list) = 
        if List.length dimensions <> 5 then failwith "read5d requires 5 dimensions"
        let r = m.Content.[name] :?> csmatio.types.MLDouble
        //This is ordered in the opposite way from previous 1D arrays in this code. 
        let data = Array.init (totalElements dimensions 1) (fun i -> r.Get(i))
        
        let dataResult = Array.init dimensions.[0] (
                            fun i -> (Array.init dimensions.[1] (
                                fun j -> (Array.init dimensions.[2] (
                                    fun k -> (Array.init dimensions.[3] (
                                        fun l -> (Array.init dimensions.[4] (
                                            fun m -> data.[i+dimensions.[0]*j+dimensions.[1]*dimensions.[0]*k+dimensions.[2]*dimensions.[1]*dimensions.[0]*l+dimensions.[3]*dimensions.[2]*dimensions.[1]*dimensions.[0]*m] ) ) ) ) ) ) ) ) )

        dataResult
        
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
                    Types.supraBasalFit=false
                    Types.matlabReplicate=true
                    }

    let cloneSizeP = read5d "PScanPP"  [(Array.length lambdaRange);(Array.length rhoRange);(Array.length rRange);(Array.length t);maxN]
    let survivalP  = read4d "PSurvScanPP" [(Array.length lambdaRange);(Array.length rhoRange);(Array.length rRange);(Array.length t)]

    let result = {  Types.cloneSizeMatrix = [|cloneSizeP|]
                    Types.survivalMatrix =  [|survivalP|]
                    Types.oneDimSizeMatrix = [||] //We don't bother calculating these at present
                    Types.oneDimSurvMatrix = [||]
                    }
    //Return a complete parametersearch result
    {  input with results = Some(result) }
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

let quickOutput filename result = 
    use file = new System.IO.StreamWriter(filename, false)//new StreamWriter(filename, false)
    List.iter (fun timepoint -> List.iter (fun i -> file.WriteLine(sprintf "%f" i)) timepoint; file.WriteLine(sprintf "")) result
    file.Close()