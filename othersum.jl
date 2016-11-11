function othersum(si::Vector{Int64},s::Vector{Float64},m::Int64)
  # OTHERSUM Compute the sum of each column of the matrix, without each
  # individual entry.  This corresponds to a matrix where each entry is the
  # sum of each column but with each entry subtracted.
  # rowsum = accumarray(si,s,[m,1])
  rowsum = zeros(Float64,m)
  for i = 1:length(si)
    rowsum[si[i]] += s[i]
  end
  os = rowsum[si] - s
  return os
end