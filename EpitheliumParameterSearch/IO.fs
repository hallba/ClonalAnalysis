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
        let rec core d acc =
            match d with
            | [] -> acc
            | topDimension::rest -> core rest acc

        let r = m.Content.[name] :?> csmatio.types.MLDouble

        //First- read all of the values into a 1D array
        let rec totalElements d acc = 
            match d with 
            | [] -> acc
            | topD::rest -> totalElements rest (acc*topD)
        
        Array.init (totalElements dimensions 1) (fun i -> r.Get(i))

    let rRange = read1d "rRange"
    let lambdaRange = read1d "lambdaRange"
    let rhoRange = read1d "rhoRange" |> Array.map (fun i -> 1. - i) //Rho correction due to gamma inversion in matlab code
    let t = read1d "t"
    let timePoints = read1d "timePoints"
    let nRange = read1d "nRange"

    //could also read maxN and nBadvalues

    let cloneSizeP = readNd "PScanPP" ([lambdaRange;rhoRange;rRange;t;nRange] |> List.map (fun i -> Array.length i))
    let survivalP  = readNd "PSurvScanPP" ([lambdaRange;rhoRange;rRange;t] |> List.map (fun i -> Array.length i))

    

    rRange