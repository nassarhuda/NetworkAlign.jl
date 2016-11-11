using MatrixNetworks
include("othermaxplus.jl")
include("othersum.jl")
include("round_messages.jl")
include("column_maxmatchsum.jl")

"""
NETALIGNMR
----------
    solve the network alignment problem with Klau's algorithm
"""

function netalignmr(S::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},
                       a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                        gamma::Float64,stepm::Int64,rtype::Int64,maxiter::Int64,verbose::Bool)


  m = maximum(li)
  n = maximum(lj)
  Matching_Setup = bipartite_matching_setup(w,li,lj,m,n)
  tripi = Matching_Setup.tripi
  matm = Matching_Setup.m
  matn = Matching_Setup.n
  rp = Matching_Setup.rp
  ci = Matching_Setup.ci
  mperm = tripi[tripi.>0]

  U = spzeros(Float64,size(S,1),size(S,2))
  xbest = zeros(Float64,length(w))

  flower = 0.0                      #best lower bound on the solution
  fupper = Inf                      #best upper bound on the solution
  next_reduction_iteration = stepm  #reduction step

  hist = zeros(Float64,maxiter,7)
  iter = 1
  @inbounds for iter = 1:maxiter
    Q = (b/2)*S + U-U'
    Qt = Q'
    Qp = Qt.colptr
    Qr = Qt.rowval
    Qv = Qt.nzval
    M = size(Qt,1)
    N = size(Qt,2)
    nedges = length(li)
    
    all_matching = column_maxmatchsum(M,N,Qp,Qr,Qv,m,n,nedges,li,lj)
    q = all_matching.q
    qj = all_matching.mi
    qi = all_matching.mj
    medges = all_matching.medges
    
    SM = sparse(qi[1:medges],qj[1:medges],1,size(Q,1),size(Q,2))

    x = a*w + q
    ai = zeros(Float64,length(tripi))

    ai[tripi.>0] = x[mperm]

    M_output = MatrixNetworks.bipartite_matching_primal_dual(rp,ci,ai,matm,matn)
    mi = MatrixNetworks.edge_indicator(M_output,li,lj)
    val = M_output.weight

    # compute stats
    matchval = dot(mi,w)
    overlap = dot(mi,(S*mi)/2)
    card = M_output.cardinality
    f = a*matchval + b*overlap

    if val < fupper
        fupper = val
        next_reduction_iteration = iter+stepm
    end

    if f > flower
        flower = f
        itermark = "*"
        xbest = convert(Vector{Float64},mi)
    else
        itermark = " "
    end


    if rtype == 1
        # no work
    elseif rtype==2
        mw = S*x
        mw = a*w + b/2*mw

        ai = zeros(Float64,length(tripi))
        ai[tripi.>0] = mw[mperm]

        M_output2 = MatrixNetworks.bipartite_matching_primal_dual(rp,ci,ai,matm,matn)
        mx = MatrixNetworks.edge_indicator(M_output2,li,lj)
        card = M_output2.cardinality
        matchval = dot(w,mx)
        overlap = dot(mx,(S*mx)/2)

        f = a*matchval + b*overlap

        if f > flower
            flower = f
            itermark = "**"
            mi = mx
            xbest = mw
        end
    end

    # report on current iter
    hist[iter,1:end] = [norm(nonzeros(U),1), flower, fupper, f, matchval, card, overlap]

    
    # the below if statement causes type instability due to @printf being type instable
    if verbose
        @printf("%5s   %4i   %8.1e   %7.2f %7.2f %7.2f  %7.2f %7.2f %7i %7i\n",
            itermark, iter, norm(nonzeros(U),1),
            flower, fupper, val,
            f, matchval, card, overlap)
    end

    if iter == next_reduction_iteration
        gamma = gamma*0.5
        
        # the below if statement causes type instability due to @printf being type instable
        if verbose
            @printf("%5s   %4s   reducing step to %f\n", "", "", gamma);
        end
        if gamma < 1e-24
            break
        end
        next_reduction_iteration = iter+stepm
    end

    if (fupper-flower) < 1e-2
        break
    end

    Wtemp = sparse(1:length(mi),1:length(mi),gamma*mi)
    U = U - Wtemp*triu(SM) + tril(SM)'*Wtemp
    Utemp1 = spones(U)
    Utemp1 *= 0.5
    U = min(U,Utemp1)
    Utemp1 *= -1
    U = max(U,Utemp1)
  end
  status = zeros(Float64,2)
  
  st = (fupper-flower) < 1e-2
  status[1] = flower
  status[2] = fupper
  hist = hist[1:iter,:]
  return (xbest,st,status,hist)

end

############################
### Additional functions ###
############################

function netalignmr(S::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},
                       a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                        gamma::Float64,stepm::Int64,rtype::Int64,maxiter::Int64)
    return netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,false)
end

function netalignmr(S::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},
                       a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                        gamma::Float64,stepm::Int64,rtype::Int64)
    return netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,1000,false)
end

function netalignmr(S::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},
                       a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                        gamma::Float64,stepm::Int64)
    return netalignmr(S,w,a,b,li,lj,gamma,stepm,1,1000,false)
end

function netalignmr(S::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},
                       a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                        gamma::Float64)
    return netalignmr(S,w,a,b,li,lj,gamma,25,1,1000,false)
end

function netalignmr(S::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},
                       a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64})
    return netalignmr(S,w,a,b,li,lj,0.4,25,1,1000,false)
end
