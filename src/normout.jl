function normout(A::SparseMatrixCSC{T,Int64}) where T
  # NORMOUT Normalize the outdegrees of the matrix A.
  #
  # P = normout(A)
  #
  #   P has the same non-zero structure as A, but is normalized such that the
  #   sum of each row is 1, assuming that A has non-negative entries.
  #

  # compute the row-sums/degrees
  #   d = sum(A,2)
  #   d = squeeze(d',1)
  #   id = 1./d
  #   P = spdiagm(id)*A
  #   return P

    S = vec(sum(A,dims=2))
    ei,ej,ev = findnz(A)
    m,n = size(A)
    P = sparse(ei,ej,ev./S[ei],m,n)
    return P

end
