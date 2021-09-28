# % NETALIGN_SETUP Setup the data for the netalign codes from the problem data

function netalign_setup(A::SparseMatrixCSC{T,Int64},B::SparseMatrixCSC{T,Int64},
                        L::SparseMatrixCSC{Q,Int64}) where {T,Q}

undirected = issymmetric(A) && issymmetric(B)
Se,Le,LeWeights = make_squares(A,B,L,undirected)

li = Le[:,1]
lj = Le[:,2]
Se1 = Se[:,1]
Se2 = Se[:,2]
S = sparse(Se1,Se2,1,nnz(L),nnz(L))

return (S,LeWeights,li,lj)

end
