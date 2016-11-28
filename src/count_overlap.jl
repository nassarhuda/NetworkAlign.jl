function count_overlap(A::SparseMatrixCSC{Int64,Int64},B::SparseMatrixCSC{Int64,Int64},mi::Array{Int64,2})

A = A[li[mi],li[mi]]
B = B[lj[mi],lj[mi]]
overlap = div(nnz(A.*B),2)
return overlap

end