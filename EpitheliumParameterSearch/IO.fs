module IO

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Data.Matlab

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
    
    let matrixToArray1D (m:Matrix<float>) = 
        let l = Matrix.columnCount m
        Array.init l (fun i -> m.[0,i])

    let matrixToArray4D (m:Matrix<float>) = 
        ()

    let m = (MatlabReader.ReadAll<float>(filename,"rhoRange","lambdaRange","rRange","PScanPP","PSurvScanPP","t","tExp","timePoints") ) // :?> System.Collections.Generic.Dictionary< string, Matrix<float> > 
    
    let rRange = m.["rRange"] |> matrixToArray1D
    let lambdaRange = m.["lambdaRange"] |> matrixToArray1D
    let rhoRange = m.["rhoRange"] |> matrixToArray1D
    let t = m.["t"] |> matrixToArray1D
    let tExp = m.["tExp"] |> matrixToArray1D
    let timePoints = m.["timePoints"] |> matrixToArray1D
    
    let l = MatlabReader.List(filename)

    //This will never work because MathNet matrices are 2D :(


    m