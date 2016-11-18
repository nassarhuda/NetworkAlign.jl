using MatrixNetworks

function load_netalign_problem(probname)

if probname == "lcsh2wiki-small"
    mainloc = joinpath(Pkg.dir("NetworkAlign"),"data/lcsh2wiki-small")
    # mainloc = "../data/lcsh2wiki-small"
    location = join([mainloc,"_A.smat"])
    A = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_B.smat"])
    B = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_L.smat"])
    (rows,header) = readdlm(location;header=true)
    L = sparse(
               convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1, 
               convert(Array{Int64,1},rows[1:parse(Int,header[3]),2])+1, 
               rows[1:parse(Int,header[3]),3],
               parse(Int,header[1]), 
               parse(Int,header[2])
               )
    # L = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_S.smat"])
    S = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_li.smat"])
    (rows,header) = readdlm(location;header=true)
    li = zeros(Int64,parse(Int,header[1]))
    ri = convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1
    li[ri] = convert(Array{Int64,1},rows[1:parse(Int,header[3]),3])

    # li = MatrixNetworks.readSMAT(location)
    # li = full(vec(li))
    
    location = join([mainloc,"_lj.smat"])
    (rows,header) = readdlm(location;header=true)
    lj = zeros(Int64,parse(Int,header[1]))
    ri = convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1
    lj[ri] = convert(Array{Int64,1},rows[1:parse(Int,header[3]),3])
    
    # lj = MatrixNetworks.readSMAT(location)
    # lj = full(vec(lj))
    
    location = join([mainloc,"_lw.smat"])
    (rows,header) = readdlm(location;header=true)
    w = zeros(Float64,parse(Int,header[1]))
    ri = convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1
    w[ri] = rows[1:parse(Int,header[3]),3]
    
    # w = MatrixNetworks.readSMAT(location)
    # w = full(vec(w))

    return (S,w,li,lj,A,B,L)
    
elseif probname == "example-overlap"
    
    mainloc = joinpath(Pkg.dir("NetworkAlign"),"data/example-overlap")
    # mainloc = "../data/example-overlap"
    location = join([mainloc,"_A.smat"])
    A = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_B.smat"])
    B = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_L.smat"])
    (rows,header) = readdlm(location;header=true)
    L = sparse(
               convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1, 
               convert(Array{Int64,1},rows[1:parse(Int,header[3]),2])+1, 
               rows[1:parse(Int,header[3]),3],
               parse(Int,header[1]), 
               parse(Int,header[2])
               )
    # L = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_S.smat"])
    S = MatrixNetworks.readSMAT(location)
    
    location = join([mainloc,"_li.smat"])
    (rows,header) = readdlm(location;header=true)
    li = zeros(Int64,parse(Int,header[1]))
    ri = convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1
    li[ri] = convert(Array{Int64,1},rows[1:parse(Int,header[3]),3])

    # li = MatrixNetworks.readSMAT(location)
    # li = full(vec(li))
    
    location = join([mainloc,"_lj.smat"])
    (rows,header) = readdlm(location;header=true)
    lj = zeros(Int64,parse(Int,header[1]))
    ri = convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1
    lj[ri] = convert(Array{Int64,1},rows[1:parse(Int,header[3]),3])
    
    # lj = MatrixNetworks.readSMAT(location)
    # lj = full(vec(lj))
    
    location = join([mainloc,"_lw.smat"])
    (rows,header) = readdlm(location;header=true)
    w = zeros(Float64,parse(Int,header[1]))
    ri = convert(Array{Int64,1},rows[1:parse(Int,header[3]),1])+1
    w[ri] = rows[1:parse(Int,header[3]),3]
    
    # w = MatrixNetworks.readSMAT(location)
    # w = full(vec(w))

    return (S,w,li,lj,A,B,L)
    
else
    error("NetAlign: problem name not recognized")
end
end