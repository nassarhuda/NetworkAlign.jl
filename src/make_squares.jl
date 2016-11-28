# % MAKE_SQUARES Returns a list of all the squares between A, B, and L

function make_squares(A::SparseMatrixCSC{Int64,Int64},B::SparseMatrixCSC{Int64,Int64},
                        L::SparseMatrixCSC{Int64,Int64},undirected::Bool)


m = size(A,1)
n = size(B,1)

if undirected
    rpA = A.colptr
    ciA = A.rowval

    rpB = B.colptr
    ciB = B.rowval
else
    A = A'
    B = B'
    
    rpA = A.colptr
    ciA = A.rowval

    rpB = B.colptr
    ciB = B.rowval
end

L = L'
rpAB = L.colptr
ciAB = L.rowval
vAB = L.nzval

Se1 = Int64[]
Se2 = Int64[]

wv = zeros(Int64,n)
sqi = 0
@inbounds for i=1:m
    # label everything in to i in B
    # label all the nodes in B that are possible matches to the current i
    possible_matches = rpAB[i]:rpAB[i+1]-1
    
    # get the exact node ids via ciA[possible_matches]
    wv[ciAB[possible_matches]] = possible_matches
    
    # for each node in A, follow its neighbors in L, then in B
    @inbounds for ri1 = rpA[i]:rpA[i+1]-1
        # get the actual node index
        ip = ciA[ri1]
        if i == ip
            continue
        end
        # for node index ip, check the nodes that its related to in L
        @inbounds for ri2 = rpAB[ip]:rpAB[ip+1]-1
            jp = ciAB[ri2]
            # for each of the nodes found in L, find their exact index in B and see if they form a square
            @inbounds for ri3 = rpB[jp]:rpB[jp+1]-1
                j = ciB[ri3]
                if j == jp
                    continue
                end
                if wv[j]>0
                    # we have a square!
                    push!(Se1,ri2)
                    push!(Se2,wv[j])
                end
            end
        end
    end
    # remove labels for things in in adjacent to B
    wv[ciAB[possible_matches]] = 0
end

Le = zeros(Int64,nnz(L),3)    
LeWeights = zeros(Float64,nnz(L))    
@inbounds for i = 1:m
    j = rpAB[i]:rpAB[i+1]-1
    Le[j,1] = i
    Le[j,2] = ciAB[j]
    LeWeights[j] = vAB[j]
end

Se = zeros(Int64,length(Se1),2)
Se[:,1] = Se1
Se[:,2] = Se2
return (Se,Le,LeWeights)
end

function make_squares(A::SparseMatrixCSC{Int64,Int64},B::SparseMatrixCSC{Int64,Int64},
                        L::SparseMatrixCSC{Float64,Int64},undirected::Bool)


m = size(A,1)
n = size(B,1)

if undirected
    rpA = A.colptr
    ciA = A.rowval

    rpB = B.colptr
    ciB = B.rowval
else
    A = A'
    B = B'
    
    rpA = A.colptr
    ciA = A.rowval

    rpB = B.colptr
    ciB = B.rowval
end

L = L'
rpAB = L.colptr
ciAB = L.rowval
vAB = L.nzval

Se1 = Int64[]
Se2 = Int64[]

wv = zeros(Int64,n)
sqi = 0
@inbounds for i=1:m
    # label everything in to i in B
    # label all the nodes in B that are possible matches to the current i
    possible_matches = rpAB[i]:rpAB[i+1]-1
    
    # get the exact node ids via ciA[possible_matches]
    wv[ciAB[possible_matches]] = possible_matches
    
    # for each node in A, follow its neighbors in L, then in B
    @inbounds for ri1 = rpA[i]:rpA[i+1]-1
        # get the actual node index
        ip = ciA[ri1]
        if i == ip
            continue
        end
        # for node index ip, check the nodes that its related to in L
        @inbounds for ri2 = rpAB[ip]:rpAB[ip+1]-1
            jp = ciAB[ri2]
            # for each of the nodes found in L, find their exact index in B and see if they form a square
            @inbounds for ri3 = rpB[jp]:rpB[jp+1]-1
                j = ciB[ri3]
                if j == jp
                    continue
                end
                if wv[j]>0
                    # we have a square!
                    push!(Se1,ri2)
                    push!(Se2,wv[j])
                end
            end
        end
    end
    # remove labels for things in in adjacent to B
    wv[ciAB[possible_matches]] = 0
end

Le = zeros(Int64,nnz(L),2)    
LeWeights = zeros(Float64,nnz(L))    
@inbounds for i = 1:m
    j = rpAB[i]:rpAB[i+1]-1
    Le[j,1] = i
    Le[j,2] = ciAB[j]
    LeWeights[j] = vAB[j]
end

Se = zeros(Int64,length(Se1),2)
Se[:,1] = Se1
Se[:,2] = Se2
return (Se,Le,LeWeights)
end

function build_Si(Se::Array{Int64,2},Le::Array{Int64,2})
# to build Si
Se1 = Se[:,1]
Se2 = Se[:,2]
Si = zeros(Int64,4,length(Se1))
for i = 1:length(Se1)
    Si[1,i] = Le[Se2[i],1]
    Si[2,i] = Le[Se1[i],1]
    Si[3,i] = Le[Se2[i],2]
    Si[4,i] = Le[Se1[i],2]
end
return Si
end