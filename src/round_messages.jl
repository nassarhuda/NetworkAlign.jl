function round_messages{T,F,R}(messages::Vector{Float64},
                        S::SparseMatrixCSC{T,Int64},w::Vector{F},
                        alpha::R,beta::Int64,
                        rp::Vector{Int64},ci::Vector{Int64},
                        tripi::Vector{Int64},n::Int64,
                        m::Int64,perm::Vector{Int64},
                        li::Vector{Int64},lj::Vector{Int64})
  ai = zeros(Float64,length(tripi))
  ai[tripi.>0] = messages[perm]
  M_output = MatrixNetworks.bipartite_matching_primal_dual(rp,ci,ai,m,n)
  mi = MatrixNetworks.edge_indicator(M_output,li,lj)
  matchweight = sum(w[mi.>0])#M_output.weight
  cardinality = sum(mi) #M_output.cardinality
  overlap = dot(mi,(S*mi)/2)
  f = alpha*matchweight + beta*overlap
  info = [f matchweight cardinality overlap]
  return info
end