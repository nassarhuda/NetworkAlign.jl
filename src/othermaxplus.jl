function othermaxplus(dim::Int64,li::Vector{Int64},lj::Vector{Int64},
                        lw::Vector{T},m::Int64,n::Int64) where T
  # OTHERMAXPLUS Apply the other-max-plus operator to a sparse matrix
  #
  # The other-max-plus operator applies the max-plus aggegration function
  # over the rows or columns of the matrix, where, for
  # each non-zeros, the non-zero itself cannot be the maximum.  Consequently,
  # Hence, it's the max-plus of all the "other" elements in the row or column
  # (which is controlled by dim).
  #
  # omp = othermaxplus(dim,li,lj,lw,m,n) applies the other-max-plus operator
  # to the matrix sparse(li, lj, lw, m, n).  The return value omp gives
  # the other-max-plus value for each non-zero element, e.g. the matrix
  # is sparse(li, lj, omp, m, n).

  if dim == 1
    # max-plus over cols
    i1 = lj
    i2 = li
    N = n
  else
    # max-plus over rows
    i1 = li
    i2 = lj
    N = m
  end
  # the output of the other-max-plus is either the maximum element
  # in the row or column (if the element itself isn't the maximum) or
  # the second largest element in the row or column (if the element itself
  # IS the maximum).

  dimmax1 = zeros(Float64,N)      # largest value
  dimmax2 = zeros(Float64,N)     # second largest value
  #this is correct because of the definition of the max-plus function.
  dimmaxind = zeros(Float64,N)   # index of largest value
  nedges = length(li)
  for i = 1:nedges
    if lw[i] > dimmax2[i1[i]]
      if lw[i] > dimmax1[i1[i]]
        dimmax2[i1[i]] = dimmax1[i1[i]]
        dimmax1[i1[i]] = lw[i]
        dimmaxind[i1[i]] = i2[i]
      else
        dimmax2[i1[i]] = lw[i]
      end
    end
  end
  omp = zeros(Float64,size(lw))
  for i = 1:nedges
    if i2[i] == dimmaxind[i1[i]]
      omp[i] = dimmax2[i1[i]]
    else
      omp[i] = dimmax1[i1[i]]
    end
  end
  return omp
end